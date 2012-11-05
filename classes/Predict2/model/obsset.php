<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @category Models
 * @abstract Bearing to satellite from observer
 */
class Predict2_Model_ObsSet
{
    public $az         = 0.0;  /*!< Azimuth [deg] */
    public $el         = 0.0;  /*!< Elevation [deg] */
    public $range      = 0.0;  /*!< Range [km] */
    public $range_rate = 0.0;  /*!< Velocity [km/sec] */
}
