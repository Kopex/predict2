<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @category Models
 * @abstract Geodetic position data structure.
 *
 */
class Predict2_Model_Geodetic{
  /**
  * Lattitude, rad
  *   
  * @var float
  */
  public $lat; 
  /**
  * Longitude, rad
  * 
  * @var float
  */
  public $lon;
  /**
  * Altitude, km
  *  
  * @var float
  */
  public $alt; 
  
  public $theta;
}
