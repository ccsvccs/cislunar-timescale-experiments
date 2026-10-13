Priori to March 1, 2013, SPICE binary kernel files (dexxx.bsp) by default 
included the position of the Mars planet center (499) with respect to 
the Mars barycenter (4) with a value of zero.
Since the masses of Phobos and Deimos are small, the actual offset of the
Mars planet center to the Mars berycenter is small, of order a few centimeters,
so using zero as an approximation is reasonable for most purposes.

Starting with DE430 released March 2013, the default is to not include
the location of the Mars planet center with respect to the Mars barycenter.
Only the position of the Mars system barycenter will be included.
The Mars planet center should come for a Mars satellite ephemeris file instead.
However, for the use of users for whom using a satellite ephemris file is cumbersome,
a modified version of the planetary ephemeris file will be included with
the location of the Mars planet center taken from a current satellite ephemeris file.

For DE430, a special version, de430_plus_MarsPC.bsp, includes the Mars planet center
location taken from the satellite ephemeris file mar097.

For DE440, similar to DE430, a special version, de440_plus_MarsPC.bsp, includes the 
Mars planet center location taken from the satellite ephemeris file mar097.

_____________________________________________________________

Starting with de430, a numerically integrated difference TT-TDB between
terrestrial time, TT, and barycentric dynamical time, TDB, is being made
availble on the SPICE binary kernals. This may be accessed by requesting
the state (position and velocity) of body 1000000001 with respect to
center 1000000000, with the value of TT-TDB in seconds returned as the
'x' component, and the rate of change, d(TT-TDB)/d(TDB) in seconds/second
returned as the 'Vx' component.

W. M. Folkner
15 January 2014

_____________________________________________________________

The original polynomial fit to TT-TDB for DE430 had numerical limits due to the
degree of the polynomial used. A revised fit with higher order polynomials
is provided in the file TTmTDB.de430.19feb2015.bsp . This gives a smoother fit,
to double precision numerical limits, although systematic uncertainties from
limited accuracy of GM of the bodies are larger than the numerical noise.

W. M. Folkner

01 November 2016

_____________________________________________________________

For de440 and de441, KBOs were added.  The following list gives the KBO IDs and names.

 N   ID    SPK_ID   NAME         ID     NAME
 1   8001  2136199  2003 UB313   136199 Eris
 2   8002  2136108  2003 EL61    136108 Haumea
 3   8003  2090377  2003 VB12     90377 Sedna
 4   8004  2136472  2005 FY9     136472 Makemake
 5   8005  2050000  2002 LM60     50000 Quaoar
 6   8006  2084522  2002 TC302    84522 -
 7   8007  2090482  2004 DW       90482 Orcus
 8   8008  2020000  2000 WR106    20000 Varuna
 9   8009  2055637  2002 UX25     55637 -
10   8010  2028978  2001 KX76     28978 Ixion
11   8011  2307261  2002 MS4     307261 -
12   8012  2174567  2003 MW12    174567 Varda
13   8013  3361580  2006 QH181        0 -
14   8014  3308265  2004 XR190        0 -
15   8015  2055565  2002 AW197    55565 -
16   8016  2145452  2005 RN43    145452 -
17   8017  2090568  2004 GV9      90568 -
18   8018  2208996  2003 AZ84    208996 -
19   8019  2225088  2007 OR10    225088 whitesnow
20   8020  2019521  1998 WH24     19521 Chaos
21   8021  2120347  2004 SB60    120347 Salacia
22   8022  2278361  2007 JJ43    278361 -
23   8023  3525142  2010 KZ39         0 -
24   8024  2230965  2004 XA192   230965 -
25   8025  2042301  2001 UR163    42301 -
26   8026  2455502  2003 UZ413        0 -
27   8027  3545742  2010 RF43         0 -
28   8028  2523639  2010 RE64         0 -
29   8029  2528381  2008 ST291        0 -
30   8030  3515022  2010 FX86         0 -

20 January 2021

_____________________________________________________________


