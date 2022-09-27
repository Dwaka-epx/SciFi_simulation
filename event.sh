#!/bin/bash 


# root 'macro.cxx(10,"my number")'
root -l -b -q 'event_wavy.C("gps_p0.5GeV_mu_water_BisY_toX")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_mu_water_BisY_toZ")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_el_water_BisY_toX")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_el_water_BisY_toZ")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_gm_water_BisY_toX")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_gm_water_BisY_toZ")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_pr_water_BisY_toX")' &
root -l -b -q 'event_wavy.C("gps_p0.5GeV_pr_water_BisY_toZ")' &
