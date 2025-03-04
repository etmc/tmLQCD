We use R to collect time and produce plots.

To collect time informations from the file `example_optput.out` type

```
    Rscript timing.R example_logfile_tmLQCD.log
```

The `timing.R` file will first collect the data and then rendering the file `Profile.Rmd`.
The time information is saved is an object mnl_list which is a list of mnl

```
    mnl:
 
        $name 
 
        $time 
 
        $list(sub_func)
```

the last entry is a list of `sub_func` which is an object of the type

```
    sub_func:
 
        $name 
 
        $time 
 
        $list(sub_func)
```
