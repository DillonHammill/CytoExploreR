# cyto_stats_compute

                    name OVAConc Treatment       alias     CD8      Va2    CD69
    1  Activation_29.fcs      50    Stim-D CD4 T Cells  169.80 26409.20  603.40
    2  Activation_30.fcs      50    Stim-D CD4 T Cells  206.10 23411.46  593.37
    3  Activation_31.fcs     500    Stim-D CD4 T Cells  260.09 21838.60 1709.09
    4  Activation_32.fcs     500    Stim-D CD4 T Cells  228.18 23924.25 1008.41
    5  Activation_33.fcs       0        NA CD4 T Cells     NaN      NaN     NaN
    6  Activation_29.fcs      50    Stim-D CD8 T Cells 9014.43 17402.61  335.73
    7  Activation_30.fcs      50    Stim-D CD8 T Cells 7957.74 13691.67  256.36
    8  Activation_31.fcs     500    Stim-D CD8 T Cells 9642.70 17075.62  372.26
    9  Activation_32.fcs     500    Stim-D CD8 T Cells 8987.55 15676.11  343.72
    10 Activation_33.fcs       0        NA CD8 T Cells     NaN      NaN     NaN
       Hoechst-405 Hoechst-430    CD44     CD4  CD11c
    1        76.68       53.45  753.50 2592.63 -45.64
    2        83.37       57.21  704.30 2650.00 -28.27
    3       112.29       77.33 1375.51 2415.09 -38.04
    4        83.66       65.89 1026.03 2499.86 -16.16
    5          NaN         NaN     NaN     NaN    NaN
    6        73.10       79.82  445.04    3.17  19.89
    7        75.81       78.16  387.51   -7.78   3.65
    8        82.17       74.24  391.61   20.28   1.12
    9        77.22       81.73  289.07   22.66 -14.41
    10         NaN         NaN     NaN     NaN    NaN

---

                    name OVAConc Treatment       alias     parent value
    1  Activation_29.fcs      50    Stim-D CD4 T Cells Live Cells    11
    2  Activation_30.fcs      50    Stim-D CD4 T Cells Live Cells    13
    3  Activation_31.fcs     500    Stim-D CD4 T Cells Live Cells    10
    4  Activation_32.fcs     500    Stim-D CD4 T Cells Live Cells    10
    5  Activation_33.fcs       0        NA CD4 T Cells Live Cells     0
    6  Activation_29.fcs      50    Stim-D CD8 T Cells Live Cells    23
    7  Activation_30.fcs      50    Stim-D CD8 T Cells Live Cells    25
    8  Activation_31.fcs     500    Stim-D CD8 T Cells Live Cells    21
    9  Activation_32.fcs     500    Stim-D CD8 T Cells Live Cells    19
    10 Activation_33.fcs       0        NA CD8 T Cells Live Cells     0
    11 Activation_29.fcs      50    Stim-D CD4 T Cells    T Cells    28
    12 Activation_30.fcs      50    Stim-D CD4 T Cells    T Cells    32
    13 Activation_31.fcs     500    Stim-D CD4 T Cells    T Cells    29
    14 Activation_32.fcs     500    Stim-D CD4 T Cells    T Cells    31
    15 Activation_33.fcs       0        NA CD4 T Cells    T Cells   NaN
    16 Activation_29.fcs      50    Stim-D CD8 T Cells    T Cells    60
    17 Activation_30.fcs      50    Stim-D CD8 T Cells    T Cells    62
    18 Activation_31.fcs     500    Stim-D CD8 T Cells    T Cells    60
    19 Activation_32.fcs     500    Stim-D CD8 T Cells    T Cells    58
    20 Activation_33.fcs       0        NA CD8 T Cells    T Cells   NaN

---

                   name       alias    FSC-A             PE-A
    1 Activation_29.fcs CD8 T Cells 65628.86 16144.7282059355
    2 Activation_30.fcs CD8 T Cells 64911.52 13486.6401414788
    3 Activation_31.fcs CD8 T Cells  69868.3 16889.1528808521
    4 Activation_32.fcs CD8 T Cells 69266.73 15433.7367957411
    5 Activation_33.fcs CD8 T Cells      NaN              NaN

---

    # A tibble: 5 x 6
      name              OVAConc Treatment alias       `FSC-A`      Va2
      <chr>             <chr>   <chr>     <chr>         <dbl>    <dbl>
    1 Activation_29.fcs 50      Stim-D    CD4 T Cells  57430.    1564.
    2 Activation_30.fcs 50      Stim-D    CD4 T Cells  62468. 1285928.
    3 Activation_31.fcs 500     Stim-D    CD4 T Cells  60474. 9660731.
    4 Activation_32.fcs 500     Stim-D    CD4 T Cells  57749.  871938.
    5 Activation_33.fcs 0       NA        CD4 T Cells     NA       NA 

---

                    name       alias  FUN           channel    value
    1  Activation_29.fcs CD4 T Cells mean Alexa Fluor 700-A  2635.63
    2  Activation_30.fcs CD4 T Cells    N Alexa Fluor 700-A       75
    3  Activation_31.fcs CD4 T Cells mean Alexa Fluor 700-A   2968.1
    4  Activation_32.fcs CD4 T Cells    N Alexa Fluor 700-A       90
    5  Activation_33.fcs CD4 T Cells mean Alexa Fluor 700-A  2835.83
    6  Activation_29.fcs CD4 T Cells    N Alexa Fluor 700-A       76
    7  Activation_30.fcs CD4 T Cells mean Alexa Fluor 700-A   2667.7
    8  Activation_31.fcs CD4 T Cells    N Alexa Fluor 700-A       72
    9  Activation_32.fcs CD4 T Cells mean Alexa Fluor 700-A      NaN
    10 Activation_33.fcs CD4 T Cells    N Alexa Fluor 700-A        0
    11 Activation_29.fcs CD4 T Cells mean              PE-A 26952.32
    12 Activation_30.fcs CD4 T Cells    N              PE-A       75
    13 Activation_31.fcs CD4 T Cells mean              PE-A 23172.01
    14 Activation_32.fcs CD4 T Cells    N              PE-A       90
    15 Activation_33.fcs CD4 T Cells mean              PE-A 23356.01
    16 Activation_29.fcs CD4 T Cells    N              PE-A       76
    17 Activation_30.fcs CD4 T Cells mean              PE-A 22493.69
    18 Activation_31.fcs CD4 T Cells    N              PE-A       72
    19 Activation_32.fcs CD4 T Cells mean              PE-A      NaN
    20 Activation_33.fcs CD4 T Cells    N              PE-A        0

