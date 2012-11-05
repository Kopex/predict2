<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 */

/*
 * Unit Solar
 *           Author:  Dr TS Kelso
 * Original Version:  1990 Jul 29
 * Current Revision:  1999 Nov 27
 *          Version:  1.30
 *        Copyright:  1990-1999, All Rights Reserved
 *
 *   Ported to C by: Neoklis Kyriazis  April 1 2001
 */
class Predict2_Solar extends Predict2_SGPObs
{
    Static $ra,$dec,$lat,$lon;
    /* Calculates solar position vector */
    public static function Calculate_Solar_Position($time, Predict2_Model_Vector $solar_vector)
    {
        $mjd = $time - 2415020.0;
        $year = 1900 + $mjd / 365.25;
        $T = ($mjd + Predict2_Time::Delta_ET($year) / Predict2::secday) / 36525.0;
        $M = Predict2_Math::Radians(Predict2_Math::Modulus(358.47583 + Predict2_Math::Modulus(35999.04975 * $T, 360.0)
             - (0.000150 + 0.0000033 * $T) * ($T * $T), 360.0));
        $L = Predict2_Math::Radians(Predict2_Math::Modulus(279.69668 + Predict2_Math::Modulus(36000.76892 * $T, 360.0)
             + 0.0003025 * ($T * $T), 360.0));
        $e = 0.01675104 - (0.0000418 + 0.000000126 * $T) * $T;
        $C = Predict2_Math::Radians((1.919460 - (0.004789 + 0.000014 * $T) * $T) * sin($M)
             + (0.020094 - 0.000100 * $T) * sin(2 * $M) + 0.000293 * sin(3 * $M));
        $O = Predict2_Math::Radians(Predict2_Math::Modulus(259.18 - 1934.142 * $T, 360.0));
        $Lsa = Predict2_Math::Modulus($L + $C - Predict2_Math::Radians(0.00569 - 0.00479 * sin($O)), Predict2::twopi);
        $nu = Predict2_Math::Modulus($M + $C, Predict2::twopi);
        $R = 1.0000002 * (1 - ($e * $e)) / (1 + $e * cos($nu));
        $eps = Predict2_Math::Radians(23.452294 - (0.0130125 + (0.00000164 - 0.000000503 * $T) * $T) * $T + 0.00256 * cos($O));
        $R = Predict2::AU * $R;

        $solar_vector->x = $R * cos($Lsa);
        $solar_vector->y = $R * sin($Lsa) * cos($eps);
        $solar_vector->z = $R * sin($Lsa) * sin($eps);
        $solar_vector->w = $R;
    }

    /* Calculates stellite's eclipse status and depth */
    public static function Sat_Eclipsed(Predict2_Model_Vector $pos, Predict2_Model_Vector $sol, &$depth)
    {
        $Rho   = new Predict2_Model_Vector();
        $earth = new Predict2_Model_Vector();

        /* Determine partial eclipse */
        $sd_earth = Predict2_Math::ArcSin(Predict2::xkmper / $pos->w);
        Predict2_Math::Vec_Sub($sol, $pos, $Rho);
        $sd_sun = Predict2_Math::ArcSin(Predict2::__sr__ / $Rho->w);
        Predict2_Math::Scalar_Multiply(-1, $pos, $earth);
        $delta = Predict2_Math::Angle($sol, $earth);
        $depth = $sd_earth - $sd_sun - $delta;

        if ($sd_earth < $sd_sun) {
            return 0;
        } else if ($depth >= 0) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * Finds the current location of the sun based on the observer location
     *
     * @param Predict2_Model_QTH $qth    The observer location
     * @param int         $daynum The daynum or null to use the current daynum
     *
     * @return Predict2_Model_ObsSet
     */
    public static function FindSun(Predict2_Model_QTH $qth, $daynum = null)
    {
        if ($daynum === null) {
            $daynum = Predict2_Time::get_current_daynum();
        }

        $obs_geodetic = new Predict2_Model_Geodetic();
        $obs_geodetic->lon   = $qth->lon * predict2::de2ra;
        $obs_geodetic->lat   = $qth->lat * predict2::de2ra;
        $obs_geodetic->alt   = $qth->alt / 1000.0;
        $obs_geodetic->theta = 0;

        $solar_vector = new Predict2_Model_Vector();
        $zero_vector  = new Predict2_Model_Vector();
        $solar_set    = new Predict2_Model_ObsSet();

        self::Calculate_Solar_Position($daynum, $solar_vector);
        Predict2_SGPObs::Calculate_Obs(
            $daynum,
            $solar_vector,
            $zero_vector,
            $obs_geodetic,
            $solar_set
        );

        $solar_set->az = Predict2_Math::Degrees($solar_set->az);
        $solar_set->el = Predict2_Math::Degrees($solar_set->el);
        //+++
        $solar_latlonalt = new Predict2_Model_Geodetic();
        Predict2_SGPObs::Calculate_LatLonAlt($daynum, $solar_vector, $solar_latlonalt);

        self::$lat=Predict2_Math::Degrees($solar_latlonalt->lat);
        self::$lon=360.0-Predict2_Math::Degrees($solar_latlonalt->lon);
        
        $solar_rad = new Predict2_Model_Vector();
        Predict2_SGPObs::Calculate_RADec($daynum, $solar_vector, $zero_vector, $obs_geodetic, $solar_rad);

        self::$ra=Predict2_Math::Degrees($solar_rad->x);
        self::$dec=Predict2_Math::Degrees($solar_rad->y);        
        //+++
        return $solar_set;
    }
}
