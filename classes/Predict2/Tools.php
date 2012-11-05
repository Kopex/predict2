<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @abstract    The main Predict2 class. 
 */

/*
    A limited PHP port of the gpredict program, done by Bill Shupp.
    Original notes and author information is below.  GPL2 license.
    ===============================================================



    Gpredict: Real-time satellite tracking and orbit prediction program

    Copyright (C)  2001-2009  Alexandru Csete, OZ9AEC.
    Parts are Copyright John A. Magliacane, KD2BD 1991-2003 (indicated below)

    Authors: Alexandru Csete <oz9aec@gmail.com>
             John A. Magliacane, KD2BD.

    Comments, questions and bugreports should be submitted via
    http://sourceforge.net/projects/gpredict/
    More details can be found at the project home page:

            http://gpredict.oz9aec.net/

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, visit http://www.fsf.org/
*/

class Predict2_Tools
{
    /* visibility constants */
    const SAT_VIS_NONE     = 0;
    const SAT_VIS_VISIBLE  = 1;
    const SAT_VIS_DAYLIGHT = 2;
    const SAT_VIS_ECLIPSED = 3;

    /* preferences */
    public $minEle     = 10; // Minimum elevation
    public $timeRes    = 10; // Pass details: time resolution
    public $numEntries = 20; // Pass details: number of entries
    public $threshold  = -6; // Twilight threshold

    /**
     *  Predict2 the next pass.
     *
     * This function simply wraps the get_pass function using the current time
     * as parameter.
     *
     *       Note: the data in sat will be corrupt (future) and must be refreshed
     *       by the caller, if the caller will need it later on (eg. if the caller
     *       is GtkSatList).
     *
     * @param Predict2_Sat $sat   The satellite data.
     * @param Predict2_Model_QTH $qth   The observer data.
     * @param int         $maxdt The maximum number of days to look ahead.
     *
     * @return Predict2_Model_Pass Pointer instance or NULL if no pass can be found.
     */
    public function get_next_pass(Predict2_Sat $sat, Predict2_Model_QTH $qth, $maxdt)
    {
        /* get the current time and call the get_pass function */
        $now = Predict2_Time::get_current_daynum();

        return $this->get_pass($sat, $qth, $now, $maxdt);
    }

    /** 
     * Predict2 first pass after a certain time.
     * 
     * This function will find the first upcoming pass with AOS no earlier than
     * t = start and no later than t = (start+maxdt).
     *
     *       For no time limit use maxdt = 0.0
     *       the data in sat will be corrupt (future) and must be refreshed
     *       by the caller, if the caller will need it later on
     *
     * @param Predict2_Sat $sat   The satellite data.
     * @param Predict2_Model_QTH $qth   The observer's location data.
     * @param float       $start Starting time.
     * @param int         $maxdt The maximum number of days to look ahead (0 for no limit).
     * 
     * @return Predict2_Model_Pass or NULL if there was an error.
     */
    public function get_pass(Predict2_Sat $sat_in, Predict2_Model_QTH $qth, $start, $maxdt)
    {
        $aos = 0.0;    /* time of AOS */
        $tca = 0.0;    /* time of TCA */
        $los = 0.0;    /* time of LOS */
        $dt = 0.0;     /* time diff */
        $step = 0.0;   /* time step */
        $t0 = $start;
        $tres = 0.0;   /* required time resolution */
        $max_el = 0.0; /* maximum elevation */
        $pass = null;
        $detail = null;
        $done = false;
        $iter = 0;      /* number of iterations */
        /* FIXME: watchdog */

        /*copy sat_in to a working structure*/
        $sat         = clone $sat_in;
        $sat_working = clone $sat_in;

        /* get time resolution; sat-cfg stores it in seconds */
        $tres = $this->timeRes / 86400.0;

        /* loop until we find a pass with elevation > SAT_CFG_INT_PRED_MIN_EL
            or we run out of time
            FIXME: we should have a safety break
        */
        while (!$done) {
            /* Find los of next pass or of current pass */
            $los = $this->find_los($sat, $qth, $t0, $maxdt); // See if a pass is ongoing
            $aos = $this->find_aos($sat, $qth, $t0, $maxdt);
            /* sat_log_log(SAT_LOG_LEVEL_MSG, "%s:%s:%d: found aos %f and los %f for t0=%f", */
            /*          __FILE__,  */
            /*          __FUNCTION__, */
            /*          __LINE__, */
            /*          aos, */
            /*          los,  */
            /*          t0); */
            if ($aos > $los) {
                // los is from an currently happening pass, find previous aos
                $aos = $this->find_prev_aos($sat, $qth, $t0);
            }

            /* aos = 0.0 means no aos */
            if ($aos == 0.0) {
                $done = true;
            } else if (($maxdt > 0.0) && ($aos > ($start + $maxdt)) ) {
                /* check whether we are within time limits;
                    maxdt = 0 mean no time limit.
                */
                $done = true;
            } else {
                //los = find_los (sat, qth, aos + 0.001, maxdt); // +1.5 min later
                $dt = $los - $aos;

                /* get time step, which will give us the max number of entries */
                $step = $dt / $this->numEntries;

                /* but if this is smaller than the required resolution
                    we go with the resolution
                */
                if ($step < $tres) {
                    $step = $tres;
                }

                /* create a pass_t entry; FIXME: g_try_new in 2.8 */
                $pass = new Predict2_Model_Pass();

                $pass->aos      = $aos;
                $pass->los      = $los;
                $pass->max_el   = 0.0;
                $pass->aos_az   = 0.0;
                $pass->los_az   = 0.0;
                $pass->maxel_az = 0.0;
                $pass->vis      = '---';
                $pass->satname  = $sat->nickname;
                $pass->details  = array();

                /* iterate over each time step */
                for ($t = $pass->aos; $t <= $pass->los; $t += $step) {

                    /* calculate satellite data */
                    $this->predict_calc($sat, $qth, $t);

                    /* in the first iter we want to store
                        pass->aos_az
                    */
                    if ($t == $pass->aos) {
                        $pass->aos_az = $sat->az;
                        $pass->orbit  = $sat->orbit;
                    }

                    /* append details to sat->details */
                    $detail             = new Predict2_Model_PassDetail();
                    $detail->time       = $t;
                    $detail->pos->x     = $sat->pos->x;
                    $detail->pos->y     = $sat->pos->y;
                    $detail->pos->z     = $sat->pos->z;
                    $detail->pos->w     = $sat->pos->w;
                    $detail->vel->x     = $sat->vel->x;
                    $detail->vel->y     = $sat->vel->y;
                    $detail->vel->z     = $sat->vel->z;
                    $detail->vel->w     = $sat->vel->w;
                    $detail->velo       = $sat->velo;
                    $detail->az         = $sat->az;
                    $detail->el         = $sat->el;
                    $detail->range      = $sat->range;
                    $detail->range_rate = $sat->range_rate;
                    $detail->lat        = $sat->ssplat;
                    $detail->lon        = $sat->ssplon;
                    $detail->alt        = $sat->alt;
                    $detail->ma         = $sat->ma;
                    $detail->phase      = $sat->phase;
                    $detail->footprint  = $sat->footprint;
                    $detail->orbit      = $sat->orbit;
                    $detail->vis        = $this->get_sat_vis($sat, $qth, $t);

                    /* also store visibility "bit" */
                    switch ($detail->vis) {
                        case self::SAT_VIS_VISIBLE:
                            $pass->vis[0] = 'V';
                            break;
                        case self::SAT_VIS_DAYLIGHT:
                            $pass->vis[1] = 'D';
                            break;
                        case self::SAT_VIS_ECLIPSED:
                            $pass->vis[2] = 'E';
                            break;
                        default:
                            break;
                    }

                    // Using an array, no need to prepend and reverse the list
                    // as gpredict does
                    $pass->details[] = $detail;

                    // Look up apparent magnitude if this is a visible pass
                    if ($detail->vis === self::SAT_VIS_VISIBLE) {
                        $apmag = $sat->calculateApparentMagnitude($t, $qth);
                        if ($pass->max_apparent_magnitude === null || $apmag < $pass->max_apparent_magnitude) {
                            $pass->max_apparent_magnitude = $apmag;
                        }
                    }

                    /* store elevation if greater than the
                        previously stored one
                    */
                    if ($sat->el > $max_el) {
                        $max_el         = $sat->el;
                        $tca            = $t;
                        $pass->maxel_az = $sat->az;
                    }

                    /*     g_print ("TIME: %f\tAZ: %f\tEL: %f (MAX: %f)\n", */
                    /*           t, sat->az, sat->el, max_el); */
                }

                /* calculate satellite data */
                $this->predict_calc($sat, $qth, $pass->los);
                /* store los_az, max_el and tca */
                $pass->los_az = $sat->az;
                $pass->max_el = $max_el;
                $pass->tca    = $tca;

                /* check whether this pass is good */
                if ($max_el >= $this->minEle) {
                    $done = true;
                } else {
                    $done = false;
                    $t0 = $los + 0.014; // +20 min
                    $pass = null;
                }

                $iter++;
            }
        }

        return $pass;
    }

    /**
     * Calculate satellite visibility.
     *
     * @param Predict2_Sat $sat     The satellite structure.
     * @param Predict2_Model_QTH $qth     The QTH
     * @param float       $jul_utc The time at which the visibility should be calculated.
     *
     * @return int The visiblity constant, 0, 1, 2, or 3 (see above)
     */
    public function get_sat_vis(Predict2_Sat $sat, Predict2_Model_QTH $qth, $jul_utc)
    {
        /* gboolean sat_sun_status;
        gdouble  sun_el;
        gdouble  threshold;
        gdouble  eclipse_depth;
        sat_vis_t vis = SAT_VIS_NONE; */

        $eclipse_depth  = 0.0;
        $zero_vector    = new Predict2_Model_Vector();
        $obs_geodetic   = new Predict2_Model_Geodetic();

        /* Solar ECI position vector  */
        $solar_vector = new Predict2_Model_Vector();

        /* Solar observed az and el vector  */
        $solar_set = new Predict2_Model_ObsSet();

        /* FIXME: could be passed as parameter */
        $obs_geodetic->lon   = $qth->lon * predict2::de2ra;
        $obs_geodetic->lat   = $qth->lat * predict2::de2ra;
        $obs_geodetic->alt   = $qth->alt / 1000.0;
        $obs_geodetic->theta = 0;

        Predict2_Solar::Calculate_Solar_Position($jul_utc, $solar_vector);
        Predict2_SGPObs::Calculate_Obs($jul_utc, $solar_vector, $zero_vector, $obs_geodetic, $solar_set);

        if (Predict2_Solar::Sat_Eclipsed($sat->pos, $solar_vector, $eclipse_depth)) {
            /* satellite is eclipsed */
            $sat_sun_status = false;
        } else {
            /* satellite in sunlight => may be visible */
            $sat_sun_status = true;
        }

        if ($sat_sun_status) {
            $sun_el = Predict2_Math::Degrees($solar_set->el);

            if ($sun_el <= $this->threshold && $sat->el >= 0.0) {
                $vis = self::SAT_VIS_VISIBLE;
            } else {
                $vis = self::SAT_VIS_DAYLIGHT;
            }
        } else {
            $vis = self::SAT_VIS_ECLIPSED;
        }

        return $vis;
    }

    /** Find the AOS time of the next pass.
     * @author Alexandru Csete, OZ9AEC
     * @author John A. Magliacane, KD2BD
     * @param Predict2_Sat $sat   The satellite data.
     * @param Predict2_Model_QTH $qth   The observer's location (QTH) data.
     * @param float       $start The julian date where calculation should start.
     * @param int         $maxdt The upper time limit in days (0.0 = no limit)
     * @return float  The julain date of the next AOS or 0.0 if the satellite has no AOS.
     *
     * This function finds the time of AOS for the first coming pass taking place
     * no earlier that start.
     * If the satellite is currently within range, the function first calls
     * find_los to get the next LOS time. Then the calculations are done using
     * the new start time.
     *
     */
    public function find_aos(Predict2_Sat $sat, Predict2_Model_QTH $qth, $start, $maxdt)
    {
        $t = $start;
        $aostime = 0.0;


        /* make sure current sat values are
            in sync with the time
        */
        $this->predict_calc($sat, $qth, $start);

        /* check whether satellite has aos */
        if (($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_GEO) ||
            ($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_DECAYED) ||
            !$this->has_aos($sat, $qth)) {

            return 0.0;
        }

        if ($sat->el > 0.0) {
            $t = $this->find_los($sat, $qth, $start, $maxdt) + 0.014; // +20 min
        }

        /* invalid time (potentially returned by find_los) */
        if ($t < 0.1) {
            return 0.0;
        }

        /* update satellite data */
        $this->predict_calc($sat, $qth, $t);

        /* use upper time limit */
        if ($maxdt > 0.0) {

            /* coarse time steps */
            while (($sat->el < -1.0) && ($t <= ($start + $maxdt))) {
                $t -= 0.00035 * ($sat->el * (($sat->alt / 8400.0) + 0.46) - 2.0);
                $this->predict_calc($sat, $qth, $t);
            }

            /* fine steps */
            while (($aostime == 0.0) && ($t <= ($start + $maxdt))) {

                if (abs($sat->el) < 0.005) {
                    $aostime = $t;
                } else {
                    $t -= $sat->el * sqrt($sat->alt) / 530000.0;
                    $this->predict_calc($sat, $qth, $t);
                }
            }
        } else {
            /* don't use upper time limit */

            /* coarse time steps */
            while ($sat->el < -1.0) {

                $t -= 0.00035 * ($sat->el * (($sat->alt / 8400.0) + 0.46) - 2.0);
                $this->predict_calc($sat, $qth, $t);
            }

            /* fine steps */
            while ($aostime == 0.0) {

                if (abs($sat->el) < 0.005) {
                    $aostime = $t;
                } else {
                    $t -= $sat->el * sqrt($sat->alt) / 530000.0;
                    $this->predict_calc($sat, $qth, $t);
                }

            }
        }

        return $aostime;
    }

    /** 
     * SGP4SDP4 driver for doing AOS/LOS calculations.
     * @param Predict2_Sat $sat The satellite data.
     * @param Predict2_Model_QTH $qth The QTH observer location data.
     * @param float       $t   The time for calculation (Julian Date)
     */
    public function predict_calc(Predict2_Sat $sat, Predict2_Model_QTH $qth = null, $t)
    {
        $obs_set      = new Predict2_Model_ObsSet();
        $sat_geodetic = new Predict2_Model_Geodetic();
        $obs_geodetic = new Predict2_Model_Geodetic();
        $sat->jul_utc = $t;
        $sat->tsince = ($sat->jul_utc - $sat->jul_epoch) * predict2::xmnpda;

        /* call the norad routines according to the deep-space flag */
        $sgpsdp = Predict2_SGPSDP::getInstance($sat);
        if ($sat->flags & Predict2_SGPSDP::DEEP_SPACE_EPHEM_FLAG) {
            $sgpsdp->SDP4($sat, $sat->tsince);
        } else {
            $sgpsdp->SGP4($sat, $sat->tsince);
        }

        Predict2_Math::Convert_Sat_State($sat->pos, $sat->vel);

        /* get the velocity of the satellite */
        $sat->vel->w = sqrt($sat->vel->x * $sat->vel->x + $sat->vel->y * $sat->vel->y + $sat->vel->z * $sat->vel->z);
        $sat->velo = $sat->vel->w;
        
        Predict2_SGPObs::Calculate_LatLonAlt($sat->jul_utc, $sat->pos, $sat_geodetic);

        while ($sat_geodetic->lon < -predict2::pi) {
            $sat_geodetic->lon += predict2::twopi;
        }

        while ($sat_geodetic->lon > (predict2::pi)) {
            $sat_geodetic->lon -= predict2::twopi;
        }
        
        if(isset($qth)){
          $obs_geodetic->lon   = $qth->lon * predict2::de2ra;
          $obs_geodetic->lat   = $qth->lat * predict2::de2ra;
          $obs_geodetic->alt   = $qth->alt / 1000.0;
          $obs_geodetic->theta = 0;
          Predict2_SGPObs::Calculate_Obs($sat->jul_utc, $sat->pos, $sat->vel, $obs_geodetic, $obs_set);
          $sat->az = Predict2_Math::Degrees($obs_set->az);
          $sat->el = Predict2_Math::Degrees($obs_set->el);
          $sat->range = $obs_set->range;
          $sat->range_rate = $obs_set->range_rate;
        }
        $sat->ssplat = Predict2_Math::Degrees($sat_geodetic->lat);
        $sat->ssplon = Predict2_Math::Degrees($sat_geodetic->lon);
        $sat->alt = $sat_geodetic->alt;
        $sat->ma = Predict2_Math::Degrees($sat->phase);
        $sat->ma *= 256.0 / 360.0;
        $sat->phase = Predict2_Math::Degrees($sat->phase);

        /* same formulas, but the one from predict is nicer */
        //sat->footprint = 2.0 * xkmper * acos (xkmper/sat->pos.w);
        $sat->footprint = 12756.33 * acos(predict2::xkmper / (predict2::xkmper + $sat->alt));
        $age = $sat->jul_utc - $sat->jul_epoch;
        $sat->orbit = floor(($sat->tle->xno * predict2::xmnpda / predict2::twopi +
                        $age * $sat->tle->bstar * predict2::ae) * $age +
                        $sat->tle->xmo / predict2::twopi) + $sat->tle->revnum - 1;
        $radec = new Predict2_Model_Vector();
        $obs_geodetic->lon = $sat->ssplon;
        $obs_geodetic->lat = $sat->ssplat;
        $obs_geodetic->alt = -Predict2::xkmper; // earth radii
        
        Predict2_SGPObs::Calculate_RADec($sat->jul_utc, $sat->pos, $sat->vel, $obs_geodetic, $radec);
        
        $sat->ra = $radec->x;
        $sat->dec = $radec->y;
    }

    /** Find the LOS time of the next pass.
     * @author Alexandru Csete, OZ9AEC
     * @author John A. Magliacane, KD2BD
     * @param Predict2_Sat $sat The satellite data.
     * @param Predict2_Model_QTH $qth The QTH observer location data.
     * @param float       $start The time where calculation should start. (Julian Date)
     * @param int         $maxdt The upper time limit in days (0.0 = no limit)
     * @return float  The time (julian date) of the next LOS or 0.0 if the satellite has no LOS.
     *
     * This function finds the time of LOS for the first coming pass taking place
     * no earlier that start.
     * If the satellite is currently out of range, the function first calls
     * find_aos to get the next AOS time. Then the calculations are done using
     * the new start time.
     * The function has a built-in watchdog to ensure that we don't end up in
     * lengthy loops.
     *
     */
    public function find_los(Predict2_Sat $sat, Predict2_Model_QTH $qth, $start, $maxdt)
    {
        $t = $start;
        $lostime = 0.0;


        $this->predict_calc($sat, $qth, $start);

        /* check whether satellite has aos */
        if (($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_GEO) ||
            ($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_DECAYED) ||
            !$this->has_aos ($sat, $qth)) {

            return 0.0;
        }

        if ($sat->el < 0.0) {
            $t = $this->find_aos($sat, $qth, $start, $maxdt) + 0.001; // +1.5 min
        }

        /* invalid time (potentially returned by find_aos) */
        if ($t < 0.01) {
            return 0.0;
        }

        /* update satellite data */
        $this->predict_calc($sat, $qth, $t);

        /* use upper time limit */
        if ($maxdt > 0.0) {

            /* coarse steps */
            while (($sat->el >= 1.0) && ($t <= ($start + $maxdt))) {
                $t += cos(($sat->el - 1.0) * predict2::de2ra) * sqrt($sat->alt) / 25000.0;
                $this->predict_calc($sat, $qth, $t);
            }

            /* fine steps */
            while (($lostime == 0.0) && ($t <= ($start + $maxdt)))  {

                $t += $sat->el * sqrt($sat->alt) / 502500.0;
                $this->predict_calc($sat, $qth, $t);

                if (abs($sat->el) < 0.005) {
                    $lostime = $t;
                }
            }
        } else {
        /* don't use upper limit */

            /* coarse steps */
            while ($sat->el >= 1.0) {
                $t += cos(($sat->el - 1.0) * predict2::de2ra) * sqrt($sat->alt) / 25000.0;
                $this->predict_calc($sat, $qth, $t);
            }

            /* fine steps */
            while ($lostime == 0.0) {

                $t += $sat->el * sqrt($sat->alt) / 502500.0;
                $this->predict_calc($sat, $qth, $t);

                if (abs($sat->el) < 0.005)
                    $lostime = $t;
            }
        }

        return $lostime;
    }

    /** 
     * Find AOS time of current pass.
     * @param Predict2_Sat $sat   The satellite to find AOS for.
     * @param Predict2_Model_QTH $qth   The ground station.
     * @param float       $start Start time, prefereably now.
     * @return The time of the previous AOS or 0.0 if the satellite has no AOS.
     *
     * This function can be used to find the AOS time in the past of the
     * current pass.
     */
    public function find_prev_aos(Predict2_Sat $sat, Predict2_Model_QTH $qth, $start)
    {
        $aostime = $start;

        /* make sure current sat values are
            in sync with the time
        */
        $this->predict_calc($sat, $qth, $start);

        /* check whether satellite has aos */
        if (($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_GEO) ||
            ($sat->otype == Predict2_SGPSDP::ORBIT_TYPE_DECAYED) ||
            !$this->has_aos($sat, $qth)) {

            return 0.0;
        }

        while ($sat->el >= 0.0) {
            $aostime -= 0.0005; // 0.75 min
            $this->predict_calc($sat, $qth, $aostime);
        }

        return $aostime;
    }

    /** 
     * Determine whether satellite ever reaches AOS.
     * @author John A. Magliacane, KD2BD
     * @author Alexandru Csete, OZ9AEC
     * @param Predict2_Sat $sat The satellite data.
     * @param Predict2_Model_QTH $qth The observer's location data
     * @return bool true if the satellite will reach AOS, false otherwise.
     */
    public function has_aos(Predict2_Sat $sat, Predict2_Model_QTH $qth)
    {
         $retcode = false;

         /* FIXME */
         if ($sat->meanmo == 0.0) {
              $retcode = false;
         } else {

            /* xincl is already in RAD by select_ephemeris */
            $lin = $sat->tle->xincl;
            if ($lin >= predict2::pio2) {
                $lin = predict2::pi - $lin;
            }

            $sma = 331.25 * exp(log(1440.0 / $sat->meanmo) * (2.0 / 3.0));
            $apogee = $sma * (1.0 + $sat->tle->eo) - predict2::xkmper;

            if ((acos(predict2::xkmper / ($apogee + predict2::xkmper)) + ($lin)) > abs($qth->lat * predict2::de2ra)) {
                $retcode = true;
            } else {
                $retcode = false;
            }
        }

        return $retcode;
    }

    /** 
     * Predict2 passes after a certain time.
     *
     *
     * This function calculates num upcoming passes with AOS no earlier
     * than t = start and not later that t = (start+maxdt). The function will
     *  repeatedly call get_pass until
     * the number of predicted passes is equal to num, the time has reached
     * limit or the get_pass function returns NULL.
     *
     *      For no time limit use maxdt = 0.0
     *      the data in sat will be corrupt (future) and must be refreshed
     *      by the caller, if the caller will need it later on (eg. if the caller
     *      is GtkSatList).
     *
     *      Prepending to a singly linked list is much faster than appending.
     *      Therefore, the elements are prepended whereafter the GSList is
     *      reversed
     *
     *
     * @param Predict2_Sat  $sat The satellite data
     * @param Predict2_Model_QTH  $qth The observer's location data
     * @param float $start The start julian date
     * @param int   $maxdt The max # of days to look
     * @param int   $num   The max # of passes to get
     * @return array of Predict2_Model_Pass instances if found, empty array otherwise
     */
    public function get_passes(Predict2_Sat $sat, Predict2_Model_QTH $qth, $start, $maxdt, $num = 0)
    {
        $passes = array();

        /* if no number has been specified
            set it to something big */
        if ($num == 0) {
            $num = 100;
        }

        $t = $start;

        for ($i = 0; $i < $num; $i++) {
            $pass = $this->get_pass($sat, $qth, $t, $maxdt);

            if ($pass != null) {
                $passes[] = $pass;
                $t = $pass->los + 0.014; // +20 min

                /* if maxdt > 0.0 check whether we have reached t = start+maxdt
                    if yes finish predictions
                */
                if (($maxdt > 0.0) && ($t >= ($start + $maxdt))) {
                    $i = $num;
                }
            } else {
                /* we can't get any more passes */
                $i = $num;
            }
        }

        return $passes;
    }

    /**
     * Filters out visible passes and adds the visible aos, tca, los, and
     * corresponding az and ele for each.
     *
     * @param array $passes The passes returned from get_passes()
     *
     * @author Bill Shupp
     * @return array
     */
    public function filterVisiblePasses(array $passes)
    {
        $filtered = array();

        foreach ($passes as $result) {
            // Dummy check
            if ($result->vis[0] != 'V') {
                continue;
            }

            $aos    = false;
            $aos_az = false;
            $aos    = false;
            $tca    = false;
            $los_az = false;
            $max_el = 0;

            foreach ($result->details as $detail) {
                if ($detail->vis != Predict2::SAT_VIS_VISIBLE) {
                    continue;
                }
                if ($detail->el < $this->minEle) {
                    continue;
                }

                if ($aos == false) {
                    $aos       = $detail->time;
                    $aos_az    = $detail->az;
                    $aos_el    = $detail->el;
                    $tca       = $detail->time;
                    $los       = $detail->time;
                    $los_az    = $detail->az;
                    $los_el    = $detail->el;
                    $max_el    = $detail->el;
                    $max_el_az = $detail->el;
                    continue;
                }
                $los    = $detail->time;
                $los_az = $detail->az;
                $los_el = $detail->el;

                if ($detail->el > $max_el) {
                    $tca       = $detail->time;
                    $max_el    = $detail->el;
                    $max_el_az = $detail->az;
                }
            }

            if ($aos === false) {
                // Does not reach minimum elevation, skip
                continue;
            }

            $result->visible_aos       = $aos;
            $result->visible_aos_az    = $aos_az;
            $result->visible_aos_el    = $aos_el;
            $result->visible_tca       = $tca;
            $result->visible_max_el    = $max_el;
            $result->visible_max_el_az = $max_el_az;
            $result->visible_los       = $los;
            $result->visible_los_az    = $los_az;
            $result->visible_los_el    = $los_el;

            $filtered[] = $result;
        }

        return $filtered;
    }

    /**
     * Translates aziumuth degrees to compass direction:
     *
     * N (0В°), NNE (22.5В°), NE (45В°), ENE (67.5В°), E (90В°), ESE (112.5В°),
     * SE (135В°), SSE (157.5В°), S (180В°), SSW (202.5В°), SW (225В°),
     * WSW (247.5В°), W (270В°), WNW (292.5В°), NW (315В°), NNW (337.5В°)
     *
     * @param int $az The azimuth in degrees, defaults to 0
     *
     * @return string
     */
    public function azDegreesToDirection($az = 0)
    {
        $i = floor($az / 22.5);
        $m = (22.5 * (2 * $i + 1)) / 2;
        $i = ($az >= $m) ? $i + 1 : $i;

        return trim(substr('N  NNENE ENEE  ESESE SSES  SSWSW WSWW  WNWNW NNWN  ', $i * 3, 3));
    }
}
