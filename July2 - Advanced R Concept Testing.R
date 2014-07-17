# July 2, 2014

# Goals: To test comprehension of 'Advanced R' by Hadley Wickham Especially 'Data
# Structures' Chapter Can be archived later.

# rm(list=ls())

data(mtcars)

y <- 1:10
attr(y, "my_attribute") <- "This is a vector"
attr(y, "my_attribute")
str(attributes(y))
attributes(y)

# Naming Vectors
x <- c(a = 1, b = 2, c = 3)
x
x <- c(1, 2, 3)
names(x) <- c("a", "b", "c")
x

x[2] <- "c"  # Can't use values not in the levels

c(factor("a"), factor("b"))  # You can't combine factors

# Factors
x <- factor(c("a", "b", "b", "a"))
x

class(x)
levels(x)

x[2] <- "c"

sex_char <- c("m", "m", "m")
sex_factor <- factor(sex_char, levels = c("m", "f"))

table(sex_char)

table(sex_factor)

# More on Factors
z <- read.csv(text = "value\n12\n1\n.\n9")
z
z$value
typeof(z$value)
typeof(z)
as.double(z$value)
class(z$value)
as.double(as.character(z$value))

# Smarter way to read in data
z <- read.csv(text = "value\n12\n1\n.\n9", na.strings = ".")
# Would have been useful at FDP
typeof(z$value)
class(z$value)
z$value

# Comment attributes
structure(1:5, comment = "my attribute")
structure(1:5, my_attribute = "This is a vector")

# Levels and Facotrs
f1 <- factor(letters) 
