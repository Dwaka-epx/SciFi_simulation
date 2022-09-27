#!/bin/bash 


# root 'macro.cxx(10,"my number")'
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_mu_water_BisY_toX_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_mu_water_BisY_toZ_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_el_water_BisY_toX_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_el_water_BisY_toZ_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_gm_water_BisY_toX_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_gm_water_BisY_toZ_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_pr_water_BisY_toX_1")' &
root -l -b -q 'event_dwavy.C("gps_p0.5GeV_pr_water_BisY_toZ_1")' &
