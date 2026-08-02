

# Top down first
# Then Side
# Change colors

# Set Template
mol modmaterial 0 0 AOChalky
mol modstyle 0 0 Licorice 0.2 100 100
mol modcolor 0 0 ColorID 1

# Set MIP Polymers
mol modmaterial 0 1 AOChalky
# mol modstyle 0 1 Licorice 0.1 100 100
mol modstyle 0 1 Lines 1.000000
mol modselect 0 1 resname "X.*"
mol modcolor 0 1 ColorID 6

# Surface
mol addrep 1
material add copy AOChalky
material change opacity Material23 0.350000
mol modstyle 1 1 QuickSurf 1.000000 0.500000 0.500000 3.000000
mol modmaterial 0 1 Material23
mol modmaterial 1 1 Material23
mol modcolor 1 1 ColorID 6

# mol addrep 1
# mol modselect 1 1 not resname "X.*"
# mol modcolor 1 1 ColorID 16

material change shininess AOChalky 0.000000

# Display stuff, cue density kinda handles white balanc
display cuedensity 0.190000
display farclip set 15.000000

# Set projection
display projection Orthographic
# display resize 1800 1200

# General changes
display shadows on
display ambientocclusion on
display aoambient 0.55
display aodirect 0.55
display antialias on

# Colors
color Display Background white
menu color off

# Turn off axis in corner
axes location Off
scale by 0.800000

# wait until colors change to render
display update 

# Top render

# render TachyonInternal mol_top.png display %s

# rotate x by -90

# render TachyonInternal mol_side.png display %s

# rotate x by 90
# rotate y by -90

# render TachyonInternal mol_side_y.png display %s
# exit