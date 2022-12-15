
export LD_PRELOAD=/usr/lib/libprofiler.so
CPUPROFILE="myprof.log" R -f gp-probit.R
CPUPROFILE="myprofa.log" R -f gp-probit2.R
google-pprof --text /usr/lib/R/bin/exec/R myprof.log | more
google-pprof --text /usr/lib/R/bin/exec/R myprofa.log | more
google-pprof --text --focus=emcore /usr/lib/R/bin/exec/R myprof.log | more
