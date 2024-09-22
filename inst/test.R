devtools::load_all("~/Documents/Git/ggstar/")
library(ggplot2)

ggplot(mtcars, aes(x=wt, y=mpg, fill=cyl)) +
  geom_star(starshape=35,size=4)
