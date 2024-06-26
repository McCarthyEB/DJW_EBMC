#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -31.79*x up 33.06*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {< -4.08,  -7.06, -26.22>, < 12.23,  -7.06, -26.22>, Rcell pigment {Black}}
cylinder {<-12.23,   7.06, -26.22>, <  4.08,   7.06, -26.22>, Rcell pigment {Black}}
cylinder {<-12.23,   7.06,   0.00>, <  4.08,   7.06,   0.00>, Rcell pigment {Black}}
cylinder {< -4.08,  -7.06,   0.00>, < 12.23,  -7.06,   0.00>, Rcell pigment {Black}}
cylinder {< -4.08,  -7.06, -26.22>, <-12.23,   7.06, -26.22>, Rcell pigment {Black}}
cylinder {< 12.23,  -7.06, -26.22>, <  4.08,   7.06, -26.22>, Rcell pigment {Black}}
cylinder {< 12.23,  -7.06,   0.00>, <  4.08,   7.06,   0.00>, Rcell pigment {Black}}
cylinder {< -4.08,  -7.06,   0.00>, <-12.23,   7.06,   0.00>, Rcell pigment {Black}}
cylinder {< -4.08,  -7.06, -26.22>, < -4.08,  -7.06,   0.00>, Rcell pigment {Black}}
cylinder {< 12.23,  -7.06, -26.22>, < 12.23,  -7.06,   0.00>, Rcell pigment {Black}}
cylinder {<  4.08,   7.06, -26.22>, <  4.08,   7.06,   0.00>, Rcell pigment {Black}}
cylinder {<-12.23,   7.06, -26.22>, <-12.23,   7.06,   0.00>, Rcell pigment {Black}}
atom(<-10.60,   6.12, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #0 
atom(< -2.45,  -2.35, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #1 
atom(< -4.08,   0.47, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #2 
atom(< -7.34,   6.12, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #3 
atom(<  2.45,  -5.18, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #4 
atom(<  0.82,  -2.35, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #5 
atom(< -0.82,   0.47, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #6 
atom(< -2.45,   3.29, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #7 
atom(< -4.08,   6.12, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #8 
atom(<  5.71,  -5.18, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #9 
atom(<  2.45,   0.47, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #10 
atom(<  0.82,   3.29, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #11 
atom(< -0.82,   6.12, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #12 
atom(<  8.97,  -5.18, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #13 
atom(<  7.34,  -2.35, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #14 
atom(<  5.71,   0.47, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #15 
atom(<  4.08,   3.29, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #16 
atom(<  2.45,   6.12, -25.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #17 
atom(<  0.82,  -0.47, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #18 
atom(< -0.82,  -3.29, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #19 
atom(< -0.82,   2.35, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #20 
atom(< -7.34,   2.35, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #21 
atom(< -2.45,   5.18, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #22 
atom(< -2.45,  -0.47, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #23 
atom(<  7.34,  -6.12, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #24 
atom(< -2.45,  -6.12, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #25 
atom(<  5.71,  -3.29, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #26 
atom(< -4.08,   2.35, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #27 
atom(<  4.08,  -0.47, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #28 
atom(< -8.97,   5.18, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #29 
atom(<  2.45,   2.35, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #30 
atom(< -5.71,   5.18, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #31 
atom(<  0.82,   5.18, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #32 
atom(< -5.71,  -0.47, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #33 
atom(< 10.60,  -6.12, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #34 
atom(<  4.08,  -6.12, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #35 
atom(<  8.97,  -3.29, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #36 
atom(<  0.82,  -6.12, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #37 
atom(<  7.34,  -0.47, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #38 
atom(<  2.45,  -3.29, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #39 
atom(<  5.71,   2.35, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #40 
atom(< -4.08,  -3.29, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #41 
atom(<  4.08,   5.18, -24.60>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #42 
atom(<  5.71,  -3.29, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #43 
atom(< -2.45,  -0.47, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #44 
atom(<  2.45,  -3.29, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #45 
atom(<  4.08,  -0.47, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #46 
atom(< -7.34,   2.35, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #47 
atom(<  0.82,  -6.12, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #48 
atom(<  2.45,   2.35, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #49 
atom(<  0.82,  -0.47, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #50 
atom(< -4.08,   2.35, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #51 
atom(<  0.82,   5.18, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #52 
atom(< -5.71,  -0.47, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #53 
atom(< -0.82,   2.35, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #54 
atom(< 10.60,  -6.12, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #55 
atom(< -2.45,  -6.12, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #56 
atom(< -5.71,   5.18, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #57 
atom(<  8.97,  -3.29, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #58 
atom(< -2.45,   5.18, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #59 
atom(< -0.82,  -3.29, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #60 
atom(<  7.34,  -0.47, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #61 
atom(< -8.97,   5.18, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #62 
atom(<  7.34,  -6.12, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #63 
atom(<  5.71,   2.35, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #64 
atom(<  4.08,  -6.12, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #65 
atom(< -4.08,  -3.29, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #66 
atom(<  4.08,   5.18, -22.60>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #67 
atom(<  0.82,   3.29, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #68 
atom(<  2.45,  -5.18, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #69 
atom(<-10.60,   6.12, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #70 
atom(< -4.08,   6.12, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #71 
atom(< -0.82,   6.12, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #72 
atom(< -8.97,   3.29, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #73 
atom(< -5.71,   3.29, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #74 
atom(<  0.82,  -2.35, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #75 
atom(<  8.97,  -5.18, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #76 
atom(<  5.71,  -5.18, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #77 
atom(< -2.45,  -2.35, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #78 
atom(< -5.71,  -2.35, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #79 
atom(<  7.34,  -2.35, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #80 
atom(< -4.08,  -5.18, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #81 
atom(<  4.08,  -2.35, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #82 
atom(< -0.82,   0.47, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #83 
atom(<  5.71,   0.47, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #84 
atom(< -7.34,   6.12, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #85 
atom(< -0.82,  -5.18, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #86 
atom(<  2.45,   0.47, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #87 
atom(<  4.08,   3.29, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #88 
atom(< -4.08,   0.47, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #89 
atom(< -2.45,   3.29, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #90 
atom(< -7.34,   0.47, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #91 
atom(<  2.45,   6.12, -21.98>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #92 
atom(< -0.82,   0.47, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #93 
atom(< -2.45,  -2.35, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #94 
atom(< -8.97,   3.29, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #95 
atom(<  2.45,  -5.18, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #96 
atom(<  8.97,  -5.18, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #97 
atom(<  4.08,  -2.35, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #98 
atom(< -5.71,   3.29, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #99 
atom(< -2.45,   3.29, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #100 
atom(< -0.82,  -5.18, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #101 
atom(<  7.34,  -2.35, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #102 
atom(<-10.60,   6.12, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #103 
atom(<  2.45,   0.47, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #104 
atom(< -5.71,  -2.35, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #105 
atom(<  0.82,  -2.35, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #106 
atom(<  5.71,   0.47, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #107 
atom(< -4.08,   6.12, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #108 
atom(< -4.08,   0.47, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #109 
atom(<  0.82,   3.29, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #110 
atom(< -7.34,   6.12, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #111 
atom(<  4.08,   3.29, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #112 
atom(< -7.34,   0.47, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #113 
atom(< -4.08,  -5.18, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #114 
atom(<  5.71,  -5.18, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #115 
atom(< -0.82,   6.12, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #116 
atom(<  2.45,   6.12, -19.98>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #117 
atom(< 10.60,  -6.12, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #118 
atom(< -2.45,  -6.12, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #119 
atom(< -2.45,   5.18, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #120 
atom(<  4.08,  -0.47, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #121 
atom(< -4.08,   2.35, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #122 
atom(< -7.34,   2.35, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #123 
atom(<  8.97,  -3.29, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #124 
atom(<  0.82,  -0.47, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #125 
atom(<  4.08,  -6.12, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #126 
atom(< -5.71,  -0.47, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #127 
atom(<  2.45,   2.35, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #128 
atom(<  7.34,  -6.12, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #129 
atom(<  7.34,  -0.47, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #130 
atom(< -2.45,  -0.47, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #131 
atom(< -8.97,   5.18, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #132 
atom(< -0.82,  -3.29, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #133 
atom(< -0.82,   2.35, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #134 
atom(<  0.82,   5.18, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #135 
atom(<  5.71,   2.35, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #136 
atom(< -5.71,   5.18, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #137 
atom(<  5.71,  -3.29, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #138 
atom(<  2.45,  -3.29, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #139 
atom(<  0.82,  -6.12, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #140 
atom(< -4.08,  -3.29, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #141 
atom(<  4.08,   5.18, -19.35>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #142 
atom(<  7.34,  -6.12, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #143 
atom(<  2.45,  -3.29, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #144 
atom(< -0.82,   2.35, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #145 
atom(<  8.97,  -3.29, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #146 
atom(<  2.45,   2.35, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #147 
atom(< -4.08,   2.35, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #148 
atom(<  0.82,  -6.12, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #149 
atom(< -2.45,  -0.47, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #150 
atom(< -0.82,  -3.29, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #151 
atom(<  5.71,  -3.29, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #152 
atom(<  7.34,  -0.47, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #153 
atom(<  4.08,  -6.12, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #154 
atom(<  0.82,   5.18, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #155 
atom(< -7.34,   2.35, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #156 
atom(< -2.45,   5.18, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #157 
atom(<  0.82,  -0.47, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #158 
atom(< -2.45,  -6.12, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #159 
atom(<  5.71,   2.35, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #160 
atom(< -4.08,  -3.29, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #161 
atom(<  4.08,  -0.47, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #162 
atom(< 10.60,  -6.12, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #163 
atom(< -5.71,   5.18, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #164 
atom(< -8.97,   5.18, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #165 
atom(< -5.71,  -0.47, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #166 
atom(<  4.08,   5.18, -17.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #167 
atom(<  7.34,  -2.35, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #168 
atom(<  4.08,  -2.35, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #169 
atom(< -4.08,   6.12, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #170 
atom(< -5.71,   3.29, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #171 
atom(< -4.08,   0.47, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #172 
atom(< -5.71,  -2.35, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #173 
atom(< -0.82,   6.12, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #174 
atom(<  0.82,  -2.35, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #175 
atom(<  5.71,   0.47, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #176 
atom(<-10.60,   6.12, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #177 
atom(< -8.97,   3.29, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #178 
atom(<  2.45,   0.47, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #179 
atom(< -2.45,   3.29, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #180 
atom(<  5.71,  -5.18, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #181 
atom(< -4.08,  -5.18, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #182 
atom(<  8.97,  -5.18, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #183 
atom(<  4.08,   3.29, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #184 
atom(<  2.45,  -5.18, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #185 
atom(< -0.82,  -5.18, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #186 
atom(< -7.34,   0.47, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #187 
atom(< -7.34,   6.12, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #188 
atom(<  0.82,   3.29, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #189 
atom(< -0.82,   0.47, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #190 
atom(< -2.45,  -2.35, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #191 
atom(<  2.45,   6.12, -16.73>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #192 
atom(< -5.71,   3.29, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #193 
atom(< -0.82,   6.12, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #194 
atom(<  5.71,  -5.18, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #195 
atom(<  2.45,   0.47, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #196 
atom(< -4.08,  -5.18, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #197 
atom(< -5.71,  -2.35, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #198 
atom(<  5.71,   0.47, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #199 
atom(< -0.82,   0.47, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #200 
atom(< -0.82,  -5.18, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #201 
atom(< -8.97,   3.29, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #202 
atom(<-10.60,   6.12, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #203 
atom(<  8.97,  -5.18, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #204 
atom(< -4.08,   6.12, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #205 
atom(<  0.82,  -2.35, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #206 
atom(<  0.82,   3.29, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #207 
atom(<  4.08,   3.29, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #208 
atom(<  4.08,  -2.35, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #209 
atom(< -2.45,  -2.35, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #210 
atom(< -7.34,   0.47, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #211 
atom(<  2.45,  -5.18, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #212 
atom(< -4.08,   0.47, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #213 
atom(<  7.34,  -2.35, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #214 
atom(< -2.45,   3.29, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #215 
atom(< -7.34,   6.12, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #216 
atom(<  2.45,   6.12, -14.74>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #217 
atom(<  7.34,  -0.47, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #218 
atom(< -5.71,  -0.47, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #219 
atom(< 10.60,  -6.12, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #220 
atom(< -0.82,   2.35, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #221 
atom(< -0.82,  -3.29, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #222 
atom(< -7.34,   2.35, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #223 
atom(< -4.08,  -3.29, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #224 
atom(<  0.82,  -0.47, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #225 
atom(<  5.71,   2.35, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #226 
atom(<  7.34,  -6.12, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #227 
atom(<  0.82,  -6.12, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #228 
atom(<  8.97,  -3.29, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #229 
atom(<  4.08,  -0.47, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #230 
atom(<  0.82,   5.18, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #231 
atom(<  2.45,  -3.29, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #232 
atom(< -4.08,   2.35, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #233 
atom(<  4.08,  -6.12, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #234 
atom(<  4.08,   5.18, -14.11>, 1.09, rgb <0.49, 0.50, 0.69>, 0.0, ase2) // #235 
atom(< -1.13,  -1.52, -12.01>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #236 
atom(< -1.13,  -1.52, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #237 
atom(<  0.82,  -2.65, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #238 
atom(< -1.13,   0.72, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #239 
atom(<  2.76,  -1.52, -12.01>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #240 
atom(<  2.76,  -1.52, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #241 
atom(<  0.82,  -0.40, -12.01>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #242 
atom(<  0.82,   1.85, -12.01>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #243 
atom(<  0.82,   1.85, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #244 
atom(<  2.76,   0.72, -10.51>, 1.24, rgb <0.00, 0.41, 0.52>, 0.0, ase2) // #245 
