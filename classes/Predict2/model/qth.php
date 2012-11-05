<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @category  Models
 * @abstract Параметры станци наблюдения.
 */

class Predict2_Model_QTH
{
  /**
  * @var string Name, eg. callsign
  */
  public $name;   /*!< Name, eg. callsign. */
  /**
  * @var string Location, eg City, Country.
  */
  public $loc;    /*!< Location, eg City, Country. */
  /**
  * @var string Описание.
  */
  public $desc;   /*!< Short description. */
  /**
  * @var float Latitude in dec. deg. North.
  */
  public $lat;    /*!< Latitude in dec. deg. North. */
  /**
  * @var float Longitude in dec. deg. East.
  */
  public $lon;    /*!< Longitude in dec. deg. East. */
  /**
  * @var float Высота над уровнем моря, метры.
  */
  public $alt;    /*!< Altitude above sea level in meters. */
  /**
  * @var mixed QRA locator
  */
  public $qra;    /*!< QRA locator */
  /**
  * @var string Weather station code (4 chars).
  */
  public $wx;     /*!< Weather station code (4 chars). */
  /**
  * @var mixed  Raw data from cfg file.
  */
  public $data;   /*!< Raw data from cfg file. */
}
