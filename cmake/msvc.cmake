# disable warning C4018 '<': signed/unsigned mismatch
# reason: MSVC OpenMP requires signed loop var
add_compile_options(/wd4018)

# disable warning C4700 uninitialized local variable 'c'
# reason: no need to initialise for ensure_precision(., *c);
add_compile_options(/wd4700)

# disable warning C4244	'=': conversion between b* and i*, u*
# reason: explicitly designed algorithm (rounding to integral type, e.g.)
add_compile_options(/wd4244)
