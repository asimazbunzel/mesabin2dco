TODO
---

[ ] improve donor switch: not just consider overflowing star but also which one is near
   to having R~RL

[ ] after changing white-dwarf to point-mass, apply f_mt = 1.0 but limit accretion to
   Eddington limit onto the WD

[ ] search in Hurley et al. 2002, MNRAS, 329, 897 and decide if we use circulariztion
   timescale of a convective or radiative envelope after a natal-kick

[ ] update natal-kick read function to avoid reading the first line as it will have
   a header in a future version of the natal-kicks tool

[ ] write better documentation

[ ] find a good test-case in the database

[ ] pass number of row to read from natal-kick filename info. as mesa assumes that
   each directory is unique per run, we need to trick it to believe so by setting
   the MESA_TEMP_CACHES_DIR to the number of row (maybe hidden and delete @ end of run)

[x] define the log_directory control for ce and cc data in order to avoid issues
   when both ce_data and cc_data folders do not exist. in order to do so, modify the
   `c*_ctrls_io.f90` files such that after its definition a call to mkdir is done
   (remember to mkdir from the utils module)
