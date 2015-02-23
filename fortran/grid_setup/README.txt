To use normal executable:
*edit p_save_grid.F
*compile with "make"
*run with "./grid_setup.sh"

To use with python
- easier to define new grids

Compile f2py function:
*if necessary, edit save_grid_f2py.F
*compile with "make py"
 (don't usually need to recompile if changing grid)

In python:
*Type:
import grid_setup as gst

*Decide on integer GRID_OPT (which grid)
*Decide on integer TEST
   0: save to proper place where rest of WIM2d code can access it
   1: save to test directory, where rest of WIM2d code cannot access it

*Type:
gst.grid_setup(GRID_OPT=GRID_OPT,TEST=TEST)
