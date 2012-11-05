<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 * @category Models
 * @abstract Geodetic position data structure.
 */
class Predict2_Model_DeepArg
{
    /* Used by dpinit part of Deep() */
    public $eosq;
    public $sinio;
    public $cosio;
    public $betao;
    public $aodp;
    public $theta2;
    public $sing;
    public $cosg;
    public $betao2;
    public $xmdot;
    public $omgdot;
    public $xnodot;
    public $xnodp;

    /* Used by dpsec and dpper parts of Deep() */
    public $xll;
    public $omgadf;
    public $xnode;
    public $em;
    public $xinc;
    public $xn;
    public $t;

    /* Used by thetg and Deep() */
    public $ds50;
}