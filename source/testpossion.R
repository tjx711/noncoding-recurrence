require(ggplot2)

p <- read.csv("http://www.ats.ucla.edu/stat/data/poisson_sim.csv")

p <- within(p, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", "Vocational"))
  id <- factor(id)
})

summary(p)

with(p, tapply(num_awards, prog, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))


ggplot(p, aes(num_awards, fill = prog)) +
  geom_histogram(binwidth=.5, position="dodge")

indexes <- which(mutcounts<10)
mutcounts <- mutcounts[indexes]
features <- features[indexes,]
ggplot(features, aes(mutcounts))+geom_histogram(binwidth=.5,position="dodge")