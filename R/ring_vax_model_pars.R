# for computing reactive vax efficacy model values
# this assumes the reactive vax efficacy response
#  - follows a logistic
#  - that 100% efficacy is achieved at 7 days
#  - that the mid point in development of efficacy is 3.5 days

ref <- function(x, k, mid=3.5) 1/(1+exp(-k*(x-3.5)))

ff <- function(
  x, k, mid=3.5, offset=ref(0,k,mid), stretch=1/(1-2*offset)
) return(stretch*(ref(x,k,mid)-offset))

require(data.table); require(cowplot)

xs <- data.table(expand.grid(x=seq(0,7,by=.1), k=seq(0.75,2.25,by=0.5)))

ggplot(xs) + aes(x=x, y=ff(x, k), color=factor(k)) +
  geom_line() +
  scale_x_continuous("days since vaccination") +
  scale_y_continuous("vaccine efficacy") +
  scale_color_discrete("logistic k") +
  theme_minimal()

# going w/ 1.25
offset <- ref(0, k=1.25)
stretch <- 1/(1-2*offset)
