<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @abstract    Contains satellite data and related methods.
 */

class Predict2_Sat extends Predict2_SGPObs 
{
    // Fifth root of a hundred, used for magnitude calculation
    const POGSONS_RATIO = 2.5118864315096;

    public $name     = null;
    public $nickname = null;
    public $website  = null;

    public $tle      = null;   /*!< Keplerian elements */
    public $flags    = 0;      /*!< Flags for algo ctrl */
    public $sgps     = null;
    public $dps      = null;
    public $deep_arg = null;
    public $pos      = null;   /*!< Raw position and range */
    public $vel      = null;   /*!< Raw velocity */

    /*** FIXME: REMOVE */
    public $bearing = null;   /*!< Az, El, range and vel */
    public $astro   = null;   /*!< Ra and Decl */
    /*** END */

    /* time keeping fields */
    public $jul_epoch = null;
    public $jul_utc   = null;
    public $tsince    = null;
    public $aos       = null;    /*!< Next AOS. */
    public $los       = null;    /*!< Next LOS */

    public $az         = null;   /*!< Azimuth [deg] */
    public $el         = null;   /*!< Elevation [deg] */
    public $range      = null;   /*!< Range [km] */
    public $range_rate = null;   /*!< Range Rate [km/sec] */
    public $ra         = null;   /*!< Right Ascension [deg] */
    public $dec        = null;   /*!< Declination [deg] */
    public $ssplat     = null;   /*!< SSP latitude [deg] */
    public $ssplon     = null;   /*!< SSP longitude [deg] */
    public $alt        = null;   /*!< altitude [km] */
    public $velo       = null;   /*!< velocity [km/s] */
    public $ma         = null;   /*!< mean anomaly */
    public $footprint  = null;   /*!< footprint */
    public $phase      = null;   /*!< orbit phase */
    public $meanmo     = null;   /*!< mean motion kept in rev/day */
    public $orbit      = null;   /*!< orbit number */
    public $otype      = null;   /*!< orbit type. */

    public function __construct(Predict2_Model_TLE $tle)
    {
        $headerParts    = explode(' ', $tle->header);
        $this->nickname = $headerParts[0];
        $this->name     = array_shift($headerParts)."\n".implode(' ',$headerParts);
        $this->tle      = $tle;
        $this->pos      = new Predict2_Model_Vector();
        $this->vel      = new Predict2_Model_Vector();
        $this->sgps     = new Predict2_Model_SGSDPStatic();
        $this->deep_arg = new Predict2_Model_DeepArg();
        $this->dps      = new Predict2_Model_DeepStatic();
        
        $this->select_ephemeris();
        $this->sat_data_init_sat($this);
    }

    /* Selects the apropriate ephemeris type to be used */
    /* for predictions according to the data in the TLE */
    /* It also processes values in the tle set so that  */
    /* they are apropriate for the sgp4/sdp4 routines   */
    public function select_ephemeris()
    {
        /* Preprocess tle set */
        $this->tle->xnodeo *= predict2::de2ra;
        $this->tle->omegao *= predict2::de2ra;
        $this->tle->xmo    *= predict2::de2ra;
        $this->tle->xincl  *= predict2::de2ra;
        $temp = predict2::twopi / predict2::xmnpda / predict2::xmnpda;

        /* store mean motion before conversion */
        $this->meanmo       = $this->tle->xno;
        $this->tle->xno     = $this->tle->xno * $temp * predict2::xmnpda;
        $this->tle->xndt2o *= $temp;
        $this->tle->xndd6o  = $this->tle->xndd6o * $temp / predict2::xmnpda;
        $this->tle->bstar  /= predict2::ae;

        /* Period > 225 minutes is deep space */
        $dd1 = predict2::xke / $this->tle->xno;
        $dd2 = predict2::tothrd;
        $a1 = pow($dd1, $dd2);
        $r1 = cos($this->tle->xincl);
        $dd1 = 1.0 - $this->tle->eo * $this->tle->eo;
        $temp = predict2::ck2 * 1.5 * ($r1 * $r1 * 3.0 - 1.0) / pow($dd1, 1.5);
        $del1 = $temp / ($a1 * $a1);
        $ao = $a1 * (1.0 - $del1 * (predict2::tothrd * 0.5 + $del1 *
                                 ($del1 * 1.654320987654321 + 1.0)));
        $delo = $temp / ($ao * $ao);
        $xnodp = $this->tle->xno / ($delo + 1.0);

        /* Select a deep-space/near-earth ephemeris */
        if (predict2::twopi / $xnodp / predict2::xmnpda >= .15625) {
            $this->flags |= Predict2_SGPSDP::DEEP_SPACE_EPHEM_FLAG;
        } else {
            $this->flags &= ~Predict2_SGPSDP::DEEP_SPACE_EPHEM_FLAG;
        }
    }

    /** Initialise satellite data.
     * @param sat The satellite to initialise.
     * @param qth Optional QTH info, use (0,0) if NULL.
     *
     * This function calculates the satellite data at t = 0, ie. epoch time
     * The function is called automatically by gtk_sat_data_read_sat.
     */
    public function sat_data_init_sat(Predict2_Sat $sat, Predict2_Model_QTH $qth = null)
    {
        $obs_geodetic = new Predict2_Model_Geodetic();
        $obs_set = new Predict2_Model_ObsSet();
        $sat_geodetic = new Predict2_Model_Geodetic();
        /* double jul_utc, age; */

        $jul_utc = Predict2_Time::Julian_Date_of_Epoch($sat->tle->epoch); // => tsince = 0.0
        $sat->jul_epoch = $jul_utc;

        /* initialise observer location */
        if ($qth != null) {
            $obs_geodetic->lon = $qth->lon * Predict2::de2ra;
            $obs_geodetic->lat = $qth->lat * Predict2::de2ra;
            $obs_geodetic->alt = $qth->alt / 1000.0;
            $obs_geodetic->theta = 0;
        }
        else {
            $obs_geodetic->lon = 0.0;
            $obs_geodetic->lat = 0.0;
            $obs_geodetic->alt = 0.0;
            $obs_geodetic->theta = 0;
        }

        /* execute computations */
        $sdpsgp = Predict2_SGPSDP::getInstance($sat);
        if ($sat->flags & Predict2_SGPSDP::DEEP_SPACE_EPHEM_FLAG) {
            $sdpsgp->SDP4($sat, 0.0);
        } else {
            $sdpsgp->SGP4($sat, 0.0);
        }

        /* scale position and velocity to km and km/sec */
        Predict2_Math::Convert_Sat_State($sat->pos, $sat->vel);

        /* get the velocity of the satellite */
        $sat->vel->w = sqrt($sat->vel->x * $sat->vel->x + $sat->vel->y * $sat->vel->y + $sat->vel->z * $sat->vel->z);
        $sat->velo = $sat->vel->w;
        Predict2_SGPObs::Calculate_Obs($jul_utc, $sat->pos, $sat->vel, $obs_geodetic, $obs_set);
        Predict2_SGPObs::Calculate_LatLonAlt($jul_utc, $sat->pos, $sat_geodetic);

        while ($sat_geodetic->lon < -predict2::pi) {
            $sat_geodetic->lon += predict2::twopi;
        }

        while ($sat_geodetic->lon > predict2::pi) {
            $sat_geodetic->lon -= predict2::twopi;
        }

        $sat->az = Predict2_Math::Degrees($obs_set->az);
        $sat->el = Predict2_Math::Degrees($obs_set->el);
        $sat->range = $obs_set->range;
        $sat->range_rate = $obs_set->range_rate;
        $sat->ssplat = Predict2_Math::Degrees($sat_geodetic->lat);
        $sat->ssplon = Predict2_Math::Degrees($sat_geodetic->lon);
        $sat->alt = $sat_geodetic->alt;
        $sat->ma = Predict2_Math::Degrees($sat->phase);
        $sat->ma *= 256.0 / 360.0;
        $sat->footprint = 2.0 * predict2::xkmper * acos (predict2::xkmper/$sat->pos->w);
        $age = 0.0;
        $sat->orbit = floor(($sat->tle->xno * predict2::xmnpda / predict2::twopi +
                                   $age * $sat->tle->bstar * predict2::ae) * $age +
                                  $sat->tle->xmo / predict2::twopi) + $sat->tle->revnum - 1;

        /* orbit type */
        $sat->otype = $sat->get_orbit_type($sat);
    }

    public function get_orbit_type(Predict2_Sat $sat)
    {
         $orbit = Predict2_SGPSDP::ORBIT_TYPE_UNKNOWN;

         if ($this->geostationary($sat)) {
              $orbit = Predict2_SGPSDP::ORBIT_TYPE_GEO;
         } else if ($this->decayed($sat)) {
              $orbit = Predict2_SGPSDP::ORBIT_TYPE_DECAYED;
         } else {
              $orbit = Predict2_SGPSDP::ORBIT_TYPE_UNKNOWN;
         }

         return $orbit;
    }


    /** Determinte whether satellite is in geostationary orbit.
     * @author John A. Magliacane, KD2BD
     * @param sat Pointer to satellite data.
     * @return TRUE if the satellite appears to be in geostationary orbit,
     *          FALSE otherwise.
     *
     * A satellite is in geostationary orbit if
     *
     *     fabs (sat.meanmotion - 1.0027) < 0.0002
     *
     * Note: Appearantly, the mean motion can deviate much more from 1.0027 than 0.0002
     */
    public function geostationary(Predict2_Sat $sat)
    {
         if (abs($sat->meanmo - 1.0027) < 0.0002) {
              return true;
         } else {
              return false;
        }
    }


    /** 
     * @author John A. Magliacane, KD2BD
     * @author Alexandru Csete, OZ9AEC
     * @param sat Pointer to satellite data.
     * @return TRUE if the satellite appears to have decayed, FALSE otherwise.
     * @version Modified version of the predict code but it is not tested.
     *
     * A satellite is decayed if
     *
     *    satepoch + ((16.666666 - sat.meanmo) / (10.0*fabs(sat.drag))) < "now"
     *
     */
    public function decayed(Predict2_Sat $sat)
    {
        /* tle.xndt2o/(twopi/xmnpda/xmnpda) is the value before converted the
           value matches up with the value in predict 2.2.3 */
        /*** FIXME decayed is treated as a static quantity.
             It is time dependent. Also sat->jul_utc is often zero
             when this function is called
        ***/
        if ($sat->jul_epoch + ((16.666666 - $sat->meanmo) /
                               (10.0 * abs($sat->tle->xndt2o / (predict2::twopi / predict2::xmnpda / predict2::xmnpda)))) < $sat->jul_utc) {
              return true;
        } else {
              return false;
        }
    }

    /**
     * Experimental attempt at calculating apparent magnitude.  Known intrinsic
     * magnitudes are listed inside the function for now.
     *
     * @param float       $time The daynum the satellite is calculated for
     * @param Predict2_Model_QTH $qth  The observer location
     *
     * @return null on failure, float otherwise
     */
    public function calculateApparentMagnitude($time, Predict2_Model_QTH $qth)
    {
        // Recorded intrinsic magnitudes and their respective
        // illumination and distance from heavens-above.com
        static $intrinsicMagnitudes = array(
            '25544' => array(
                'mag'      => -1.3,
                'illum'    => .5,
                'distance' => 1000,
            )
        );

        // Return null if we don't have a record of the intrinsic mag
        if (!isset($intrinsicMagnitudes[$this->tle->catnr])) {
            return null;
        }
        $imag = $intrinsicMagnitudes[$this->tle->catnr];

        // Convert the observer's geodetic info to radians and km so
        // we can compare vectors
        $observerGeo      = new Predict2_Model_Geodetic();
        $observerGeo->lat = Predict2_Math::Radians($qth->lat);
        $observerGeo->lon = Predict2_Math::Radians($qth->lon);
        $observerGeo->alt = $qth->alt * 1000;

        // Now determine the sun and observer positions
        $observerPos      = new Predict2_Model_Vector();
        $observerVel      = new Predict2_Model_Vector();
        $solarVector      = new Predict2_Model_Vector();
        Predict2_Solar::Calculate_Solar_Position($time, $solarVector);
        Predict2_SGPObs::Calculate_User_PosVel($time, $observerGeo, $observerPos, $observerVel);

        // Determine the solar phase and and thus the percent illumination
        $observerSatPos = new Predict2_Model_Vector();
        Predict2_Math::Vec_Sub($this->pos, $observerPos, $observerSatPos);
        $phaseAngle = Predict2_Math::Degrees(Predict2_Math::Angle($solarVector, $observerSatPos));
        $illum      = $phaseAngle / 180;

        $illuminationChange            = $illum / $imag['illum'];
        $inverseSquareOfDistanceChange = pow(($imag['distance'] / $this->range), 2);
        $changeInMagnitude             = log(
            $illuminationChange * $inverseSquareOfDistanceChange,
            self::POGSONS_RATIO
        );

        return $imag['mag'] - $changeInMagnitude;
    }
}
