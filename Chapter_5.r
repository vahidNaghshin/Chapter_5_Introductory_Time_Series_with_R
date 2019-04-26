set.seed(1)
z <- w <- rnorm(100, sd = 25)
for (t in 2:100) z[t] <- 0.5 * z[t - 1] + w[t]
Time <- 1:100
print(3**2)
x <- 70 + 2 * Time + 3* (Time**2)+ z

png(
  "quad1.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
plot(x, xlab = "time", type = "l")
dev.off()
browseURL("quad1.png")

 x.lm <- lm(x ~ Time+ I(Time^2))
 coef(x.lm)
 confint(x.lm)
 sqrt(diag(vcov(x.lm)))


 png(
  "quad1corrLm.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
pacf(resid(x.lm))
dev.off()
browseURL("quad1corrLm.png")

library(nlme)
x.gls <- gls(x ~ Time+ I(Time^2), cor = corAR1(0.44))
coef(x.gls)
confint(x.gls)
sqrt(diag(vcov(x.gls)))

www <-file.path(getwd(), 'global.dat')
Global <- scan(www)
Global.ts <- ts(Global, st = c(1856, 1), end = c(2005, 12), fr = 12)
temp <- window(Global.ts, start = 1970)

TIME <- (time(temp) - mean(time(temp)))/sd(time(temp))
SIN <- COS <- matrix(nr = length(TIME), nc = 6)
for (i in 1:6) {
     COS[, i] <- cos(2 * pi * i * time(temp))
     SIN[, i] <- sin(2 * pi * i * time(temp))
 }

temp.lm1 <- lm(temp ~ TIME + I(TIME^2) +
                        COS[,1] + SIN[,1] + COS[,2] + SIN[,2] +
                        COS[,3] + SIN[,3] + COS[,4] + SIN[,4] +
                        COS[,5] + SIN[,5] + COS[,6] + SIN[,6])
coef(temp.lm1)/sqrt(diag(vcov(temp.lm1)))

temp.lm2 <- lm(temp ~ TIME + SIN[, 1] + SIN[, 2])
coef(temp.lm2)/sqrt(diag(vcov(temp.lm2)))

png(
  "quad1corrtemp.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
pacf(resid(temp.lm2))
dev.off()
browseURL("quad1corrtemp.png")

temp.lm2 <- lm(temp ~ TIME + SIN[, 1] + SIN[, 2])
coef(temp.lm2)/sqrt(diag(vcov(temp.lm2)))
sqrt(diag(vcov(temp.lm2)))

library(nlme)
temp.gls <-gls(temp ~ TIME + SIN[, 1] + SIN[, 2], cor = corAR1(0.5))
coef(temp.gls)/sqrt(diag(vcov(temp.gls)))
sqrt(diag(vcov(temp.gls)))
# SIN[,2] is an explanatory variable that looks significant in temp.lm2 but not in temp.gls



data_cbe_path <-file.path(getwd(), 'cbe.dat')
CBE <- read.table(data_cbe_path, header = T)
class(CBE)
Elec.ts <- ts(CBE[, 3], start = 1958, freq = 12)

png(
  "elec.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
plot(Elec.ts, type='l')
dev.off()
browseURL("elec.png")

Seas <- cycle(Elec.ts)
Time <- time(Elec.ts)
print(Time)  
elec.lm <- lm(log(Elec.ts) ~ Time + I(Time^2) + factor(Seas))
step(elec.lm)

SIN <- COS <- matrix(nr = length(Time), nc = 6)
for (i in 1:6) {
     COS[, i] <- cos(2 * pi * i * Time)
     SIN[, i] <- sin(2 * pi * i * Time)
 }

elec.lm1 <- lm(log(Elec.ts) ~ Time + I(Time^2) +
                        COS[,1] + SIN[,1] + COS[,2] + SIN[,2] +
                        COS[,3] + SIN[,3] + COS[,4] + SIN[,4] +
                        COS[,5] + SIN[,5] + COS[,6] + SIN[,6])
step(elec.lm1)


png(
  "elec_res.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
acf(resid(elec.lm1))
dev.off()
browseURL("elec_res.png")

Elec_res.ar <- ar(resid(elec.lm1))

png(
  "elec_res_ar.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
acf(Elec_res.ar$res[-(1:Elec_res.ar$order)], main="residual error correlogram", lag = 50)
dev.off()
browseURL("elec_res_ar.png")

Time <- time(ts(start = 1991, end = c(2000, 12), fr = 12))
SIN <- COS <- matrix(nr = length(Time), nc = 6)
for (i in 1:6) {
     COS[, i] <- cos(2 * pi * i * Time)
     SIN[, i] <- sin(2 * pi * i * Time)
 }

 new.dat <- data.frame(TIME = Time, SIN = SIN,
      COS = COS)

ar.pred = predict(Elec_res.ar, n.ahead=120)
log.pred <- predict(elec.lm1, new.dat)
elec.pred = exp(log.pred + ar.pred$pred + 0.5*0.0003933)
elec.pred.ts = ts(elec.pred, st=1991, fr=12)




png(
  "elec_res_ar_pred.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
ts.plot(Elec.ts, elec.pred.ts, lty=1:2)
dev.off()
browseURL("elec_res_ar_pred.png")


www <-file.path(getwd(), 'Fontdsdt.dat')
Fontdsdt.dat <- read.table(www, header=T)
adflow.ts <- ts(Fontdsdt.dat$adflow, frequency = 12)

# print(adflow.ts)

png(
  "fontadst.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
plot(adflow.ts, type='l')
dev.off()
browseURL("fontadst.png")


Seas <- cycle(adflow.ts)
Time <- time(adflow.ts)
# print(Time)
adflow.lm <- lm(adflow.ts ~ 0 + Time + factor(Seas))
coef(adflow.lm)

png(
  "fontadst_res.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)
pacf(resid(adflow.lm))
dev.off()
browseURL("fontadst_res.png")

adflow.res.ar <- ar(resid(adflow.lm))

png(
  "ar_res_hist.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)

w <- rexp(1000)
hist(adflow.res.ar$res[-(1:adflow.res.ar$order)], prob = T)
dev.off()
browseURL("ar_res_hist.png")

adflow.pred.ar <- predict(adflow.res.ar, n.ahead = 120)

new.t <- seq(73, len = 10 * 12, by = 1/12)
new.dat <- data.frame(Time = new.t, Seas = rep(1:12, 10))
adflow.pred <- predict(adflow.lm, new.dat)

adflow.pred.ts <- ts(adflow.pred, st=73, fr=12)
adflow.pred.ar.ts <- ts(adflow.pred.ar$pred, st=73, fr=12)

png(
  "font_pred.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)

ts.plot(adflow.ts, adflow.pred.ts+adflow.pred.ar.ts, lty=1:2)
dev.off()
browseURL("font_pred.png")

adflow.lm <- lm(log(adflow.ts) ~ 0 + Time + factor(Seas))
coef(adflow.lm)

adflow.res.ar <- ar(resid(adflow.lm))

adflow.pred.ar <- predict(adflow.res.ar, n.ahead = 120)

new.t <- seq(73, len = 10 * 12, by = 1/12)
new.dat <- data.frame(Time = new.t, Seas = rep(1:12, 10))
adflow.pred <- predict(adflow.lm, new.dat)


exp.adflow.ts <- ts(exp(adflow.pred+adflow.pred.ar$pred), st=73, fr=12)

png(
  "font_pred_log.png",
  width     = 5.25,
  height    = 3.75,
  units     = "in",
  res       = 700,
  pointsize = 4
)

ts.plot(adflow.ts, exp.adflow.ts, lty=1:2)
dev.off()
browseURL("font_pred_log.png")
