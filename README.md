# Lid-Driven-Cavity

SIMPLE(Semi-Implicit Method for Pressure-Linked Equation) algorithm implemented for driven cavity problem on collocated grid in C++ and post-processed in Python. This approach solves the issue of checkerboard oscillations which was originally faced by developers. Implemention is slightly more difficult in comparison to staggered grid but the code has extensive comments which makes it easily readable. In addition to that, I have also outlined some useful resources in order to tackle the problem. It is expected that the user has experience with staggered formulation before approaching this problem. 

# Results(Re=100)


<p float="left">
  <img src="https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/u.eps" width="400" >
  <img src="https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/v.eps" width="400"> 
</p>

<p float="left">
  <img src="https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/p.eps" width="400" >
  <img src="https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/u_centerline.eps" width="400"> 
</p>
  
<image src= "https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/v_centerline.eps" width=400>
<image src= "https://github.com/deepmorzaria/Lid-Driven-Cavity/blob/main/results/streamlines.eps" width=400>




# Resources

[Computational Fluid Dynamics by Dr. Sumon Chakraborty](https://www.youtube.com/playlist?list=PL3zvA_WajfGBi-0-A9goGqB0cbe5-aU4N)

[CFD lectures by Dr. Sandip Mazumdar](https://www.youtube.com/playlist?list=PLVuuXJfoPgT4gJcBAAFPW7uMwjFKB9aqT)

[Introduction to Computational Fluid Dynamics by H K Versteeg and W Malalasekera](http://ftp.demec.ufpr.br/disciplinas/TM702/Versteeg_Malalasekera_2ed.pdf)
