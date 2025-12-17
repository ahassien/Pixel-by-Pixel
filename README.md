![Olive Trees](Olive_trees_generative_art_for_CHANCE.png)

# Pixel-by-Pixel
A generative art challenge was issued in the April 2023 issue of CHANCE. For this challenge, readers were asked to use their brains and computational skills to create generative art. This challenge encouraged contestants to showcase the skills and knowledge of data science while highlighting the creative applications of the discipline. 

My generative interpretation of Vincent Van Gogh's "Olive Trees" won the student category and was featured on the cover of CHANCE magazine Volume 37 Issue 2. 

I used flow fields to create dynamic curves that emulate the style of not only Van Gogh’s “Olive Trees” but also his other works such as “Starry Night.” In my artwork, the soothing color palettes for the sky, mountains, trees, and grass were based on true colors from the painting.

The R code uses internal functions from the aRtsy package as well as ggplot2. I used a cubic noise generator to initialize the angles of the flow field. I also used custom sine waves with random jitters to smooth out the transitions between the horizontal color palettes.

The current version of the code takes approximately 25 minutes to run, but a simplified version of the image can be produced by reducing the number in the “lines” variable.
