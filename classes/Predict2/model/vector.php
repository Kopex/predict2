<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @category  Models
 * @abstract General three-dimensional vector structure.
 */
class Predict2_Model_Vector 
{
    public $x = 0;
    public $y = 0;
    public $z = 0;
    public $w = 0;
    
    function __construct($x=0,$y=0,$z=0,$w=0){
      $this->x=(double)$x;
      $this->y=(double)$y;
      $this->z=(double)$z;
      $this->w=(double)$w;
    }
}
