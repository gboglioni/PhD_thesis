# Generate Figures in Thesis as at 06/09/2022

# Load required function
source("create_PIBD_figures.R")

# 'triangles' example: 2D
# No colors
pdf(file = "ex_triangles_black.pdf", width = 13, height = 7)
p.i.figures(type = 'triangles', n=2000, my.col = 'black', pt.size = 1.5, plot.type = '2D.black', var.ind = c(1,2,3))
dev.off()

# With colors
pdf(file = "ex_triangles_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'triangles', n=2000, pt.size = 1.5, plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

pdf(file = "ex_triangles_U1_U3.pdf", width = 8.5, height = 8)
p.i.figures(type = 'triangles', n=2000, pt.size = 1.5, plot.type = '2D.colour', var.ind = c(1,3,2))
dev.off()

# 'triangles' example: 3D
p.i.figures(type = 'triangles', my.col = 'brown3', n=2000, pt.size = 3/4, plot.type = '3D', type.3d = 's')

my.view <- structure(c(0.8, -0.21, 0.50, 
                       0, 0.5, 0.17, -0.8, 
                       0, 0.09, 0.96, 0.25, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('ex_triangles_3D.png', fmt = 'png')

# my.view <- par3d("userMatrix") # (Save the current viewpoint)
# dput(my.view) 3 (output coordinates of my.view)

# Continuous Bernstein: 2D
# With colors
pdf(file = "ex_Bernstein_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'bernstein', n=2000, pt.size = 2, plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

# Continuous Bernstein: 3D
p.i.figures(type = 'bernstein', my.col = 'brown3', n=2000, pt.size = 3/4, plot.type = '3D', type.3d = 's')
my.view <- structure(c(0.832010328769684, 0.244331777095795, -0.498057425022125, 
                       0, -0.546996653079987, 0.510985612869263, -0.663090288639069, 
                       0, 0.0924861133098602, 0.824133157730103, 0.558793962001801, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('ex_bernstein_3D.png', fmt = 'png')

# Tetrandron: 2D
# With colors
pdf(file = "ex_tetra_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'tetra', n=2000, pt.size = 2,  plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

# Tetraedron: 3D
p.i.figures(type = 'tetra', my.col = 'brown3', n=2000, pt.size = 3/4, plot.type = '3D', type.3d = 's')
my.view <- structure(c(0.191600948572159, 0.811443746089935, 0.488823115825653, 
                       0, 0.879559814929962, -0.346859931945801, 0.157996997237206, 
                       0, 0.3271704018116, 0.423253536224365, -0.840247631072998, 0, 
                       0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('ex_tetra_3D.png', fmt = 'png')

# Cosine: 2D
# With colors
pdf(file = "ex_cosine_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'cosine', n=2000, pt.size = 2, plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

pdf(file = "ex_cosine_U1_U3.pdf", width = 8.5, height = 8)
p.i.figures(type = 'cosine', n=2000, pt.size = 2, plot.type = '2D.colour', var.ind = c(1,3,2))
dev.off()

# Cosine: 3D
p.i.figures(type = 'cosine', my.col = 'brown3', n=2000, pt.size = 3/4, plot.type = '3D', type.3d = 's')
my.view.cosine <- structure(c(0.720872461795807, -0.148704543709755, 0.626362264156342, 
                              0, 0.588248670101166, 0.408877521753311, -0.636932194232941, 
                              0, -0.171804502606392, 0.894675612449646, 0.402788460254669, 
                              0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view.cosine)
rgl.snapshot('ex_cosine_3D.png', fmt = 'png')


# "copula" PIBD: 2D
pdf(file = "ex_copulaPIBD_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'copula', alpha =1, n=2000, pt.size = 2, plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

# "copula": 3D
p.i.figures(type = 'copula', my.col = 'brown3', n=2000,  pt.size = 3/4, plot.type = '3D', type.3d = 's')
my.view.copula <- structure(c(0.157347485423088, 0.821067869663239, 0.484963119029999, 
                              0, 0.714926958084106, -0.448799669742584, 0.454267263412476, 
                              0, 0.640993297100067, 0.288431525230408, -0.705786645412445, 
                              0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view.copula)
rgl.snapshot('ex_copulaPIBD_3D.png', fmt = 'png')


# "mixture" of 'triangles' with mutual independence PIBD: 2D
pdf(file = "ex_mix_triangles_U1_U2.pdf", width = 8.5, height = 8)
p.i.figures(type = 'triangles', alpha =0.9, n=2000, pt.size = 2, plot.type = '2D.colour', var.ind = c(1,2,3))
dev.off()

# "mixture" of 'triangles' with mutual independence PIBD: 2D
p.i.figures(type = 'triangles', alpha = 0.9, my.col = 'brown3', n=2000, pt.size = 3/4, plot.type = '3D', type.3d='s')

my.view <- structure(c(0.8, -0.21, 0.50, 
                       0, 0.5, 0.17, -0.8, 
                       0, 0.09, 0.96, 0.25, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('ex_mix_triangles_3D.png', fmt = 'png')


###### MATRICES of 3X3 scatterplots #######
# Triangles Example
pdf(file = "ex_mat_triangles.pdf", width = 10, height = 10)
p.i.figures(type = 'triangles', n=2000, pt.size = 0.75, plot.type='2D.matrix', var.ind = c(1,2,3), legend=F, margins = T)
dev.off()

# Bernstein Example
pdf(file = "ex_mat_bernstein.pdf", width = 10, height = 10)
p.i.figures(type = 'bernstein', n=2000, pt.size = 0.75, plot.type='2D.matrix', var.ind = c(1,2,3), legend=F, margins = T)
dev.off()


#################################################
######      3D FIGURES: EMPIRICAL FIT      ######
#################################################
# (GBB: modified 04/07/2022 for type.3d = 's', pt.size  should be around 1, versus 8 if type.3d = 'p')
# Independence, size n = 1000 
p.i.figures(type = 'indep', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D',
            emp.colors = T, scale.col = F, rad=0.176, type.3d = 'p')
my.view1 <- structure(c(-0.0723488107323647, 0.838174343109131, 0.475740522146225, 
                        0, 0.783040642738342, -0.25382798910141, 0.491242527961731, 0, 
                        0.580098807811737, 0.432969123125076, -0.684268176555634, 0, 
                        0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view1)
rgl.snapshot('empirical_indep_n1000_r01_grid.png', fmt = 'png')


# Bernstein, size n = 1000
p.i.figures(type = 'bernstein', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
par3d(userMatrix = my.view1)
rgl.snapshot('empirical_bernstein_n1000_r0176.png', fmt = 'png')

# Bernstein, size n = 1000, RELATIVE scale
p.i.figures(type = 'bernstein', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = T, rad=0.176)
par3d(userMatrix = my.view1)
rgl.snapshot('empirical_bernstein_n1000_r0176_relative.png', fmt = 'png')

# Triangles, size n = 1000
p.i.figures(type = 'triangles', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
my.view2 <- structure(c(0.8, -0.21, 0.50, 
                        0, 0.5, 0.17, -0.8, 
                        0, 0.09, 0.96, 0.25, 
                        0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view2)
rgl.snapshot('empirical_triangles_n1000_r0176.png', fmt = 'png')

# Tetraedron, size n = 1000
p.i.figures(type = 'tetra', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
my.view3 <- structure(c(0.191600948572159, 0.811443746089935, 0.488823115825653, 
                        0, 0.879559814929962, -0.346859931945801, 0.157996997237206, 
                        0, 0.3271704018116, 0.423253536224365, -0.840247631072998, 0, 
                        0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view3)
rgl.snapshot('empirical_tetra_n1000_r0176.png', fmt = 'png')


# Cosine, size n = 1000
p.i.figures(type = 'cosine', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
par3d(userMatrix = my.view.cosine)
rgl.snapshot('empirical_cosine_n1000_r0176.png', fmt = 'png')

# Copula, size n = 1000
p.i.figures(type = 'copula', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
my.view4 <- structure(c(0.198320806026459, 0.858878195285797, 0.39635443687439, 
                        0, 0.89115434885025, -0.089692085981369, -0.341613471508026, 
                        0, -0.270542323589325, 0.448185324668884, -0.847429394721985, 
                        0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view4)
rgl.snapshot('empirical_copula_n1000_r0176.png', fmt = 'png')

# Copula, size n = 1000, ABOVE only
p.i.figures(type = 'copula', n = 1000, my.col = c("whitesmoke", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8,  plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.176)
par3d(userMatrix = my.view4)
rgl.snapshot('empirical_copula_n1000_r0176_above.png', fmt = 'png')

# Copula, size n = 1000, RELATIVE scale
p.i.figures(type = 'copula', n = 1000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 8, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = T, rad=0.176)
par3d(userMatrix = my.view4)
rgl.snapshot('empirical_copula_n1000_r0176_relative.png', fmt = 'png')


###### GRID VERSIONS ######

# Grid method:indep
p.i.figures(type = 'indep', n = 5000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 1/2, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.103, grid=T, adjust=T)
par3d(userMatrix = my.view1)
rgl.snapshot('empirical_indep_n1000_r0103_grid.png', fmt = 'png')

# my.view <- par3d("userMatrix") # (Save the current viewpoint)
# dput(my.view) 3 (output coordinates of my.view)
p.i.figures(type = 'triangles', n = 5000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 1, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.103, grid = T, adjust.pts.size = T)
my.view <- structure(c(0.529750883579254, -0.21016538143158, 0.78057324886322, 
                       0, 0.780429780483246, 0.176978632807732, -0.527738511562347, 
                       0, -0.0249058678746223, 0.958719611167908, 0.269140124320984, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('empirical_triangles_n1000_r0103_grid.png', fmt = 'png')

#Grid version:bernstein
p.i.figures(type = 'bernstein', n = 5000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 1/2, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.103, grid.method=T, adjust.pts.size = T, g=8)
par3d(userMatrix = my.view1)
rgl.snapshot('empirical_bernstein_n1000_r0103_grid.png', fmt = 'png')

# Tetra
p.i.figures(type = 'tetra', n = 5000, my.col = c("blue", "whitesmoke", "red"), 
            alpha = 1, pt.size = 1/4, plot.type = '3D', type.3d = 'p',
            emp.colors = T, scale.col = F, rad=0.103, grid = T, adjust=T)
my.view <- structure(c(0.542121946811676, 0.600589692592621, 0.528673589229584, 
                       0, 0.740638673305511, -0.233158558607101, -0.562131524085999, 
                       0, -0.22662553191185, 0.744351983070374, -0.621917366981506, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('empirical_tetra_n1000_r0103_grid.png', fmt = 'png')


# 4D Visualisation
p.i.figures(type = 'bernstein4D', alpha = 2 , n=3000, pt.size = 10, my.col = c("blue", "white", "red"), plot.type='3D', type.3d = 'p')
my.view <- structure(c(0.67372727394104, 0.350477069616318, -0.597793698310852, 
                       0, -0.038545735180378, 0.819088339805603, 0.496496796607971, 
                       0, 0.710596799850464, -0.329546004533768, 0.615346431732178, 
                       0, 0, 0, 0, 1), .Dim = c(4L, 4L))
par3d(userMatrix = my.view)
rgl.snapshot('ex_bernstein4D.png', fmt = 'png')
