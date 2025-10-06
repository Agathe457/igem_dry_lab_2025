;
; Phytophthora/AMP ABM (AND = Elicitin ∧ ROS)
; Two-phase defense: baseline containment + pathogen-specific amplification
;

patches-own [
  chemical              ; amount of AMP present
  phyto                 ; phyto mass blob
  elicitin              ; pathogen-derived PAMP produced with phyto
  ros                   ; reactive oxygen species (from AMP/phyto/glucans)
  glucan                ; β-1,3/1,6-glucan fragments released upon AMP attack
  toxin                 ; kill switch toxin

]

turtles-own [
  m                     ; mRNA level of the output gene
  y                     ; protein level of the output gene - defines on/off
  on?                   ; if on then secreting amp iif y >= y - threshold
  active?               ; activity sensor when no phyto present
]

globals[
  ; circuit parameters
  time-step
  a0 amax                ; transcription leak & max transcription gain
  beta dm dp             ; translation, mRNA decay, protein decay
  K1 K2 n1 n2            ; Hill parameters: K1/n1 elicitin, K2/n2 for ROS
  y-threshold            ; ON threshold on Y

  ;secretion in two steps
  q                      ; baseline secretion of AMP, alway on if secretion ON
  q_boost                ; extra secretion when AND gate is ON during attack

  ; dehydration (from 0 to 1) which simultaneously reduces diffusion,
  ; increases evaporation, and reduces transcriptional gain (proxy for stress)
  diffusion-rate
  evaporation-rate
  diff-min diff-max
  evap-min evap-max

  ; pathogen (phyto)
  C_eff                    ; AMP level considered effective
  phyto-initial
  phyto-initial-max        ; max value when size is = 1 biomass value in blob
  phyto-kill-rate          ; AMP driven kill rate
  phyto-decay-rate         ; natural decay of pathogen
  phyto-present?           ; true if any phyto remains

  target-protected-area    ; number of patches wanted >= C_eff
  time-to-protection       ; first tick when target is met or -1 if not yet

  phyto-radius-max       ; radius of the blob
  phyto-radius
  blob-roughness           ; 0=perfect disk, >0 adds noisy edge
  rough-min rough-max

  amp-max                  ; scale halo intensity
  amp-max-attack           ; scale halo intensity in phyto presence

  ; inputs
  elicitin-prod              ; production by phyto
  elicitin-decay
  ros-prod-amp            ; ROS made in response to AMP
  ros-prod-phyto          ; ROS from pathogen stress/presence
  ros-prod-glucan         ; ROS amplified by glucan fragments
  ros-decay
  glucan-release-rate     ; fragments released when AMP attacks phyto
  glucan-decay



  chemical-initial
  chemical-add

  reshuffle-period      ; ticks between non active/active random
  next-reshuffle-tick

  ;plotting
  amp-secreted-per-tick
  mean-m mean-y

  ; toxin for kill switch
  toxin-dose         ; uniform toxin level injected on all patches
  lethal-dose        ; threshold above which turtles die instantly

]

to-report hill[ s K n ]
  let denom ((K ^ n) + (s ^ n))
  if denom = 0 [ report 0 ]
  report ((s ^ n) / denom)
end

to seed-phyto-blob
  let cx (min-pxcor + max-pxcor) / 2
  let cy (min-pycor + max-pycor) / 2
  ;; noisy disk include patch if distance < noisy radius
  ask patches [
    set phyto 0
    let d distancexy cx cy
    ;; noise multiplies the radius by a factor in [1 - blob-roughness, 1 + blob-roughness]
    let jitter (1 + blob-roughness * (random-float 2 - 1))
    if d <= phyto-radius * jitter [
      set phyto phyto-initial
    ]
  ]
end


to setup
  clear-all
  set-default-shape turtles "rod2"

  set time-step 1
  set reshuffle-period 30
  set next-reshuffle-tick reshuffle-period

  set a0 0.8
  set amax 1.2
  set beta 0.5
  set dm 0.10
  set dp 0.02
  set K1 0.5
  set K2 0.8
  set n1 2.0
  set n2 2.0
  set y-threshold 2
  set q 0.2
  set q_boost 0.8
  set C_eff 0.3
  set phyto-initial-max 1.0
  set phyto-radius-max 15
  set phyto-radius   (phyto-size * phyto-radius-max)
  set phyto-initial  (phyto-size * phyto-initial-max)
  set phyto-kill-rate 0.4
  set phyto-decay-rate 0
  set target-protected-area 200
  set chemical-initial 0.0
  set chemical-add 0.0

  ; dehydration set in slider derives diffusion and evaporation
  set diff-min 10      ; lowest diffusion when very dry
  set diff-max 95      ; highest diffusion when very wet
  set evap-min 5       ; lowest evaporation when very wet
  set evap-max 80      ; highest evaporation when very dry

  ; from current dehydration slider
  set diffusion-rate   (diff-max - (diff-max - diff-min) * dehydration)
  set evaporation-rate (evap-min + (evap-max - evap-min) * dehydration)

  set amp-max 0.5
  set amp-max-attack 1.5

  set elicitin-prod    0.02
  set elicitin-decay 0.0005
  set ros-prod-amp    0.05
  set ros-prod-phyto   0.02
  set ros-prod-glucan 0.06
  set ros-decay 0.10
  set glucan-release-rate  0.04
  set glucan-decay 0.01


  set toxin-dose      1.0
  set lethal-dose     0.5


  ask patches [
    set chemical chemical-initial
    set phyto 0
    set elicitin 0
    set ros 0
    set glucan 0
     set toxin 0
    set pcolor black
  ]

  set time-to-protection -1

  set rough-min 0.00    ; more compact when wet
  set rough-max 1    ; drying means rougher
  set blob-roughness (rough-min + (rough-max - rough-min) * dehydration)

  seed-phyto-blob

  ; bacteria is randomly scattered
  create-turtles population [
    setxy random-xcor random-ycor
    set size 2
    set m 0
    set y random-float 1
    set on? false
    set active? false
    ; initial color while pathogen exists = OFF > orange
    set color orange
  ]


  set phyto-present? any? patches with [ phyto > 0.0001 ]
  ask patches [
    ifelse phyto > 0 and phyto-initial > 0
      [ let b min list 1 (phyto / phyto-initial)
        set pcolor rgb round (120 * b) round (30 * b) round (25 * b) ]
    [ set pcolor black ]
  ]


  ; each turtle has a 50% chance of being active at start when no phyto present
  if not phyto-present? [
    ask turtles [ set active? (random-float 1 < 0.5) ]
  ]

  if plot? [
    clear-all-plots
  ]

  reset-ticks
end

;;;;;;;;;;;;;;;;;;;;;
;;; Go procedures ;;;
;;;;;;;;;;;;;;;;;;;;;

to go  ; forever button

  set phyto-present? any? patches with [ phyto > 0.0001 ]

  set amp-secreted-per-tick 0

  if not phyto-present? [
    if ticks >= next-reshuffle-tick [
      ask turtles [ set active? (random-float 1 < 0.5) ]
      set next-reshuffle-tick ticks + reshuffle-period
    ]
  ]

  ask turtles [ cell-step ]

  ; update each tick so can see changes without resetting
  set diffusion-rate   (diff-max - (diff-max - diff-min) * dehydration)
  set evaporation-rate (evap-min + (evap-max - evap-min) * dehydration)

  ; global background AMP injection
  if chemical-add > 0 [
    ask patches [ set chemical chemical + chemical-add * time-step ]
  ]

  ; dehydration reduces diffusion, increases evaporation
  ; use the dehydration-derived rates directly
  diffuse chemical (diffusion-rate / 100)

  ask patches [
    set chemical (chemical * (100 - evaporation-rate) / 100)
  ]


  ;   ; AMP-driven kill + natural decay
  ;  ask patches [
  ;    let kill phyto-kill-rate * chemical
  ;    set phyto max list 0 ( phyto - (phyto-decay-rate * phyto + kill) * time-step )
  ;  ]
  ;  ; this method not strong enough too linear, need stringer kill effect for visible results on blob


  ask patches [
    let Hkill ((chemical ^ 3) / ((0.8 ^ 3) + (chemical ^ 3) + 1e-9))  ;; MIC≈0.8, n=3
    let frac-kill 0.25 * Hkill        ;; max 25% of phyto per tick near/above MIC
    let dphyto (phyto-decay-rate * phyto + frac-kill * phyto)
    set phyto max list 0 (phyto - dphyto * time-step)
  ]

  ; input dynamics on patches
  ask patches [
    ;; elicitin from pathogen present?
    set elicitin max list 0 (elicitin
      + elicitin-prod * phyto * time-step
      - elicitin-decay * elicitin * time-step)

    ; glucan fragments released when AMP attacks phyto
    set glucan max list 0 (glucan
      + glucan-release-rate * chemical * phyto * time-step
      - glucan-decay * glucan * time-step)

    ; ROS from AMP, pathogen stress, and glucan-augmented signaling
    let amp_ros_term (ifelse-value (phyto > 1e-6) [ ros-prod-amp * chemical ] [ 0 ])
    set ros max list 0 (ros
      + (amp_ros_term + ros-prod-phyto * phyto + ros-prod-glucan * glucan) * time-step
      - ros-decay * ros * time-step)
  ]

  ; any phyto remains?
  set phyto-present? any? patches with [ phyto > 0.0001 ]


  ; plots
  let turtle-count count turtles
  let frac_on (ifelse-value (turtle-count = 0) [0]
    [ (count turtles with [ on? ]) / turtle-count ])
  let frac_active (ifelse-value (turtle-count = 0) [0]
    [ (count turtles with [ active? ]) / turtle-count ])

  let total_amp   sum [ chemical ] of patches
  let total_phyto sum [ phyto ] of patches
  let total_elic  sum [ elicitin ] of patches
  let total_ros   sum [ ros ] of patches
  let total_glu   sum [ glucan ] of patches
  let protected_area count patches with [ chemical >= C_eff ]

  set mean-m (ifelse-value (turtle-count = 0) [0] [ mean [ m ] of turtles ])
  set mean-y (ifelse-value (turtle-count = 0) [0] [ mean [ y ] of turtles ])

  ; capture first time protection target is reached
  if (time-to-protection = -1 and protected_area >= target-protected-area) [
    set time-to-protection ticks
  ]

; update plots
  if plot? [
; 1. Secretion
  set-current-plot "AMP Secretion"
  set-current-plot-pen "secreted"
  plot amp-secreted-per-tick

  ; Population states
  set-current-plot "Population states (on/off)"
  set-current-plot-pen "frac-on"      plot frac_on
  set-current-plot-pen "frac-active"  plot (ifelse-value phyto-present? [0] [frac_active])

  ; Concentrations
  set-current-plot "Field concentrations"
  set-current-plot-pen "AMP"      plot total_amp
  set-current-plot-pen "Phyto"    plot total_phyto
  set-current-plot-pen "Elicitin" plot total_elic
  set-current-plot-pen "ROS"      plot total_ros
  set-current-plot-pen "Glucan"   plot total_glu

  ]

  ;; recolor patches: show Phytophthora when present; otherwise show AMP old code does not work
  ;  ask patches [
  ;    ifelse (phyto > 0)
  ;    [ set pcolor scale-color blue phyto 0 phyto-initial ]
  ;    [ set pcolor scale-color green chemical 0.1 5 ]
  ;  ]
;
;  ask patches [
;    ;;new colouring that works
;    let g round (255 * chemical / amp-max)          ;; AMP intensity (green)
;    let b round (255 * (ifelse-value (phyto-initial > 0) [phyto / phyto-initial] [0]))    ;; Phyto intensity (blue)
;    set g max list 0 (min list 255 g)
;    set b max list 0 (min list 255 b)
;
;    let R 255 - (min list 255 (g + b))   ;; both AMP and Phyto dim red
;    let G2 255 - b                        ;; Phyto dims green
;    let B2 255 - g                        ;; AMP dims blue
;
;    set pcolor rgb R G2 B2
  ;  ]

  ask patches [
    let ampscale (ifelse-value phyto-present?
      [ amp-max-attack ]
      [ amp-max ])
    let ampN  min list 1 ((ifelse-value (ampScale > 0) [ chemical / ampScale ] [ 0 ]) ^ 0.6)
    let phyN  min list 1 (ifelse-value (phyto-initial > 0) [ phyto / phyto-initial ]     [ 0 ])
    let g round (255 * ampN)          ;; AMP contribution (green channel)
    let pR round (120 * phyN)         ;; phyto maroon R
    let pG round ( 30 * phyN)         ;; phyto maroon G
    let pB round ( 25 * phyN)         ;; phyto maroon B
    let R pR
    let G2 pG + g
    let B pB

    set pcolor rgb (max list 0 (min list 255 R))
    (max list 0 (min list 255 G))
    (max list 0 (min list 255 B))
  ]
  tick
end

to cell-step  ; turtle procedure
  ; immediate kill if local toxin >= lethal dose
  if [toxin] of patch-here >= lethal-dose [
    die
    stop
  ]

  ; if no phyto, go into radom ON/OFF = random active?
  if not phyto-present? [
    set on? active?

    if secretion? [
      ask patch-here [
        let q_idle (0.05 * q)                 ; 5% leak when OFF
        let add (ifelse-value [on?] of myself
          [ q * time-step ]               ; ON
          [ q_idle * time-step ])         ; OFF
        set chemical chemical + add
        set amp-secreted-per-tick amp-secreted-per-tick + add
      ]
    ]

    rt (random-float 6 - 3)
    if not can-move? 0.01 [ set heading towardsxy 0 0 ]
    fd 0.01
    ;; simple idle coloring
    set color (ifelse-value on? [cyan] [orange])
    stop
  ]

  ; sense inputs on current patch
  let s_elicitin [elicitin] of patch-here
  let s_ros      [ros]      of patch-here

  ; AND gate & dehydration-attenuated transcription gain
  let H1 hill s_elicitin K1 n1
  let H2 hill s_ros      K2 n2
  let gE (0.6 + (1.0 - 0.6) * (1 - dehydration))  ; moisture effect on expression gain: drier -> smaller gE, dehydration = 0 then gE=1.0

  ; ODE tau-leap with tiny uniform noise = adding noise, stochastic gene expression
  let dm_dt ((a0 * gE) + (amax * gE *  (H1 * H2))) - ((dm + 0.01) * m)
  ; F(H1,H2) = H1*H2 -> AND, if we want NAND -> change to  F(H1,H2) = 1 - (H1*H2)
  let dY_dt ((beta * m)) - ((dp + 0.01) * y)

  set m max (list 0 (m + dm_dt * time-step + (random-float 0.2 - 0.1)))
  set y max (list 0 (y + dY_dt * time-step + (random-float 0.2 - 0.1)))

 ; set on? (y >= y-threshold)   ;; ON/OFF decision OLD

  ; new
  ; probabilistic ON: quicker near the blob, slower farther away
  ifelse( y < y-threshold) [
    set on? false]
  [
    ;; p = exp(-d/6) with 0.1 floor; near center p≈1, far away p→0.1
    let d distancexy ((min-pxcor + max-pxcor) / 2) ((min-pycor + max-pycor) / 2)
    let p max list 0.1 (exp (- d / 6))
    if random-float 1 < p [ set on? true ]
  ]


;  ;; two-phase secretion: baseline containment always; boost only if AND is ON and invasion persists OLD
;  if secretion? [
;    ask patch-here [
;      let add (q * time-step)
;      if [on?] of myself and [phyto-present?] of myself [
;        set add add  + ( q_boost * time-step)
;      ]
;      set chemical chemical + add
;      set amp-secreted-per-tick amp-secreted-per-tick + add
;    ]
;  ]

  ;new
  if secretion? [
    if phyto-present? [
      ; distance-weighted secretion only during invasion
      let d distancexy ((min-pxcor + max-pxcor) / 2) ((min-pycor + max-pycor) / 2)
      let scale max list 0.15 (exp (- d / 6))   ;; stronger falloff; 0.05 floor
      ask patch-here [
        ; OFF = no baseline during invasion; ON = baseline + boost, both scaled
        let add (ifelse-value [on?] of myself
          [ (q * scale + q_boost * scale) * time-step ]
          [ 0 ])
        set chemical chemical + add
        set amp-secreted-per-tick amp-secreted-per-tick + add
      ]
    ]
  ]


  ; gentle wobble and inwards push
  rt (random-float 6 - 3) ;; uniform wobble
  if not can-move? 0.01 [ set heading towardsxy 0 0 ]
  fd 0.01

  ifelse on? [ set color orange ] [ set color cyan ]
  ; COLOURS:
  ; phyto present: ON: yellow, OFF: red
  ; phyto cleared: mission complete: cyan
  ifelse phyto-present? [
    ifelse on?
    [ set color yellow]   ; ON & pathogen present
    [ set color orange ]      ; OFF & pathogen present
  ] [
    set color cyan          ; pathogen cleared  mission complete
  ]

end

to kill-switch
  ; inject lethal toxin
  ask patches [ set toxin toxin-dose ]
end
@#$#@#$#@
GRAPHICS-WINDOW
257
10
762
516
-1
-1
7.0
1
10
1
1
1
0
0
0
1
-35
35
-35
35
1
1
1
ticks
30.0

BUTTON
37
74
117
107
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
137
74
212
107
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
34
22
224
55
population
population
0.0
200.0
184.0
1.0
1
NIL
HORIZONTAL

SLIDER
33
178
205
211
dehydration
dehydration
0
1
0.1
0.05
1
NIL
HORIZONTAL

SWITCH
809
29
912
62
plot?
plot?
0
1
-1000

SWITCH
933
30
1053
63
secretion?
secretion?
0
1
-1000

SLIDER
32
125
204
158
phyto-size
phyto-size
0
1
0.7
0.1
1
NIL
HORIZONTAL

PLOT
782
85
982
235
AMP Secretion
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"secreted" 1.0 0 -10899396 true "" ""

PLOT
1003
84
1203
234
Population states (on/off)
NIL
NIL
0.0
10.0
0.0
2.0
true
false
"" ""
PENS
"frac-on" 1.0 0 -4079321 true "" ""
"frac-active" 1.0 0 -5825686 true "" ""

PLOT
781
255
981
405
Field concentrations
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"AMP" 1.0 0 -2674135 true "" ""
"ROS" 1.0 0 -1184463 true "" ""
"Phyto" 1.0 0 -8630108 true "" ""
"Elicitin" 1.0 0 -10899396 true "" ""
"Glucan" 1.0 0 -2064490 true "" ""

BUTTON
78
240
173
273
NIL
kill-switch
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

In this project, a colony of ants forages for food. Though each ant follows a set of simple rules, the colony as a whole acts in a sophisticated way.

## HOW IT WORKS

When an ant finds a piece of food, it carries the food back to the nest, dropping a chemical as it moves. When other ants "sniff" the chemical, they follow the chemical toward the food. As more ants carry food to the nest, they reinforce the chemical trail.

## HOW TO USE IT

Click the SETUP button to set up the ant nest (in violet, at center) and three piles of food. Click the GO button to start the simulation. The chemical is shown in a green-to-white gradient.

The EVAPORATION-RATE slider controls the evaporation rate of the chemical. The DIFFUSION-RATE slider controls the diffusion rate of the chemical.

If you want to change the number of ants, move the POPULATION slider before pressing SETUP.

## THINGS TO NOTICE

The ant colony generally exploits the food source in order, starting with the food closest to the nest, and finishing with the food most distant from the nest. It is more difficult for the ants to form a stable trail to the more distant food, since the chemical trail has more time to evaporate and diffuse before being reinforced.

Once the colony finishes collecting the closest food, the chemical trail to that food naturally disappears, freeing up ants to help collect the other food sources. The more distant food sources require a larger "critical number" of ants to form a stable trail.

The consumption of the food is shown in a plot.  The line colors in the plot match the colors of the food piles.

## EXTENDING THE MODEL

Try different placements for the food sources. What happens if two food sources are equidistant from the nest? When that happens in the real world, ant colonies typically exploit one source then the other (not at the same time).

In this project, the ants use a "trick" to find their way back to the nest: they follow the "nest scent." Real ants use a variety of different approaches to find their way back to the nest. Try to implement some alternative strategies.

The ants only respond to chemical levels between 0.05 and 2.  The lower limit is used so the ants aren't infinitely sensitive.  Try removing the upper limit.  What happens?  Why?

In the `uphill-chemical` procedure, the ant "follows the gradient" of the chemical. That is, it "sniffs" in three directions, then turns in the direction where the chemical is strongest. You might want to try variants of the `uphill-chemical` procedure, changing the number and placement of "ant sniffs."

## NETLOGO FEATURES

The built-in `diffuse` primitive lets us diffuse the chemical easily without complicated code.

The primitive `patch-right-and-ahead` is used to make the ants smell in different directions without actually turning.

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Wilensky, U. (1997).  NetLogo Ants model.  http://ccl.northwestern.edu/netlogo/models/Ants.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 1997 Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

This model was created as part of the project: CONNECTED MATHEMATICS: MAKING SENSE OF COMPLEX PHENOMENA THROUGH BUILDING OBJECT-BASED PARALLEL MODELS (OBPML).  The project gratefully acknowledges the support of the National Science Foundation (Applications of Advanced Technologies Program) -- grant numbers RED #9552950 and REC #9632612.

This model was developed at the MIT Media Lab using CM StarLogo.  See Resnick, M. (1994) "Turtles, Termites and Traffic Jams: Explorations in Massively Parallel Microworlds."  Cambridge, MA: MIT Press.  Adapted to StarLogoT, 1997, as part of the Connected Mathematics Project.

This model was converted to NetLogo as part of the projects: PARTICIPATORY SIMULATIONS: NETWORK-BASED DESIGN FOR SYSTEMS LEARNING IN CLASSROOMS and/or INTEGRATED SIMULATION AND MODELING ENVIRONMENT. The project gratefully acknowledges the support of the National Science Foundation (REPP & ROLE programs) -- grant numbers REC #9814682 and REC-0126227. Converted from StarLogoT to NetLogo, 1998.

<!-- 1997 1998 MIT -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

rod
true
0
Circle -2674135 true false 90 105 90
Circle -2674135 true false 75 105 92
Circle -2674135 true false 180 150 0
Circle -2674135 true false 105 105 90
Circle -2674135 true false 102 132 66
Circle -2674135 true false 122 107 88

rod2
true
0
Circle -7500403 true true 45 105 90
Circle -7500403 true true 73 88 96
Circle -7500403 true true 78 93 85
Circle -7500403 true true 105 75 92

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
