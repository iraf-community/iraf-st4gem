# {LOAD.CL - LOAD SOME REQUIRED TASKS FROM STSDAS

if ( deftask("tables")) {
    if ( !defpac( "tables" )) {
        tables motd- 
        ;
    }
}
;

# Load some of the package so the user can access them without loading.
# Avoid printing menu, but do not change the
# default value of the menus switch.

      analysis
      fourier
      dither
      fitting
      graphics
      stplot 
      toolbox
      headers
      imgtools
      tools
      nttools

keep
