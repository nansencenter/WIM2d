1. To use normal executable:
   *compile with "make"
   *cd to place where you want to run the model
      mkdir Build; cd Build
   *copy infile_grid.txt and infile_waves.txt
   *edit infile_grid.txt and infile_waves.txt if necessary
   *Type
      "$WIM2D_PATH/grid_setup/scripts/grid_setup.sh"

2. To use with python
   - easier to define new grids (eg different LANDMASKS)

   Compile f2py function:
   *if necessary, edit save_grid_f2py.F
   *compile with "make py"
    (don't usually need to recompile if changing grid)

   *cd to place where you want to run the model
      mkdir Build; cd Build
   *Type:
      python $WIM2D_PATH/grid_setup/scripts/grid_setup.py (also uses the infile_grid.txt)
   or
      python $WIM2D_PATH/grid_setup/scripts/grid_setup.py --GRID_OPT=0
         (can also use GRID_OPT=1,2 for different hard-coded versions of the LAND_MASK and size)
