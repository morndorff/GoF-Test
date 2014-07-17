# July 2, 2014 Goals: Design 'Trapezoid' Area approximation for two pdfs

# Start: Consider two pdfs X ~ N(0,sqrt(3)) Y ~ t(3)

# Then for a sample of size 10, the approximation for the area would be:
set.seed(5)
x <- rnorm(10, 0, sqrt(3))
y <- rnorm(10, 3)

lenx <- length(x)
delta <- (1/lenx)/3
x1 <- seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx)
x <- qnorm(x1, 0, sqrt(3))
x_m_del <- x - delta
x_p_del <- x + delta
x1_m_del <- pnorm(x_m_del, 0, sqrt(3))
x1_p_del <- pnorm(x_p_del, 0, sqrt(3))
# Corners of X

y <- qt(x1, 3)
y_m_del <- y - delta
y_p_del <- y + delta
y1_m_del <- pt(y_m_del, 3)
y1_p_del <- pt(y_p_del, 3)

# Height of trapezoid is 2*delta Length of top side is abs(x1_p_del - y1_p_del)
# Length of bottom side is abs(x1_m_del - y1_m_del)
trap1 <- (abs(x1_p_del - y1_p_del) + abs(x1_m_del - y1_m_del)) * (delta)
area <- sum((abs(x1_p_del - y1_p_del) + abs(x1_m_del - y1_m_del)) * delta) 
