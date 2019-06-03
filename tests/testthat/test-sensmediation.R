context("sensmediation and more.effects")

data("RSdata")

test_that("effects are correct for mediator and outcome models without covariates", {

  medmod.nocovariates <- glm(lowered.consc ~ AF, data = RSdata, family = binomial(link = "probit"))
  outmod.nocovariates <- glm(cf.3mo ~ AF + lowered.consc, data = RSdata, family = binomial(link = "probit"))

  sensmed.nocovariates <- sensmediation(med.model = medmod.nocovariates, out.model = outmod.nocovariates, exp.name = "AF1",
                                        med.name = "lowered.consc")

  b0 <- medmod.nocovariates$coefficients["(Intercept)"]
  b1 <- medmod.nocovariates$coefficients["AF1"]
  th0 <- outmod.nocovariates$coefficients["(Intercept)"]
  th1 <- outmod.nocovariates$coefficients["AF1"]
  th2 <- outmod.nocovariates$coefficients["lowered.consc"]

  nie <- matrix((stats::pnorm(th0 + th2 + th1) - stats::pnorm(th0 + th1))*(stats::pnorm(b0 + b1)- stats::pnorm(b0)))
  nde <- matrix((stats::pnorm(th0 + th1) - stats::pnorm(th0))*(1 - stats::pnorm(b0) )  +
                  (stats::pnorm(th0 + th2 + th1) - stats::pnorm(th0 + th2))*stats::pnorm(b0))

  expect_equivalent(sensmed.nocovariates$NIE, nie)
  expect_equivalent(sensmed.nocovariates$NDE, nde)

})

test_that("effects are correct for mediator and outcome models with exactly one covariate", {

  medmod.1covariate <- glm(lowered.consc ~ AF + sex, data = RSdata, family = binomial(link = "probit"))
  outmod.1covariate <- glm(cf.3mo ~ AF + lowered.consc + sex, data = RSdata, family = binomial(link = "probit"))

  sensmed.1covariate <- sensmediation(med.model = medmod.1covariate, out.model = outmod.1covariate, exp.name = "AF1",
                                      med.name = "lowered.consc")

  b0 <- medmod.1covariate$coefficients["(Intercept)"]
  b1 <- medmod.1covariate$coefficients["AF1"]
  b2 <- medmod.1covariate$coefficients["sex1"]
  th0 <- outmod.1covariate$coefficients["(Intercept)"]
  th1 <- outmod.1covariate$coefficients["AF1"]
  th2 <- outmod.1covariate$coefficients["lowered.consc"]
  th4 <- outmod.1covariate$coefficients["sex1"]
  x <- matrix(RSdata$sex)
  x <- matrix(as.numeric(x))

  nie <- matrix(mean((stats::pnorm(th0 + th2 + th1 + x%*%th4) - stats::pnorm(th0 + th1 + x%*%th4))*(stats::pnorm(b0 + b1 + x%*%b2) - stats::pnorm(b0 + x%*%b2))))
  nde <- matrix(mean((stats::pnorm(th0 + th1 + x%*%th4) - stats::pnorm(th0 + x%*%th4))*(1 - stats::pnorm(b0 + x%*%b2) )  +
                       (stats::pnorm(th0 + th2 + th1 + x%*%th4) - stats::pnorm(th0 + th2 + x%*%th4))*stats::pnorm(b0 + x%*%b2) ))

  expect_equivalent(sensmed.1covariate$NIE, nie)
  expect_equivalent(sensmed.1covariate$NDE, nde)


})

test_that("coefficients/covariances correctly stored for mediator and outcome models with interactions involving the exposure and/or mediator", {

  medmod.interactions <- glm(lowered.consc ~ AF + sex + age.cat + AF*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))
  outmod.interactions <- glm(cf.3mo ~ sex + age.cat + AF * lowered.consc + lowered.consc*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))

  sensmed.interactions <- sensmediation(med.model = medmod.interactions, out.model = outmod.interactions, exp.name = "AF1",
                                        med.name = "lowered.consc")
  sensmed.interactions$betas
  sensmed.interactions$thetas
  sensmed.interactions$sigma.thetabeta

  glm.betas <- list(medmod.interactions$coefficients["(Intercept)"],medmod.interactions$coefficients["AF1"],
                    matrix(c(medmod.interactions$coefficients["sex1"], medmod.interactions$coefficients["age.cat70-79"],
                             medmod.interactions$coefficients["age.cat80-89"], medmod.interactions$coefficients["age.cat90-"])),
                    matrix(c(medmod.interactions$coefficients["AF1:sex1"], medmod.interactions$coefficients["AF1:age.cat70-79"],
                             medmod.interactions$coefficients["AF1:age.cat80-89"], medmod.interactions$coefficients["AF1:age.cat90-"])))

  glm.thetas <- list(outmod.interactions$coefficients["(Intercept)"],outmod.interactions$coefficients["AF1"],
                     outmod.interactions$coefficients["lowered.consc"], outmod.interactions$coefficients["AF1:lowered.consc"],
                     matrix(c(outmod.interactions$coefficients["sex1"], outmod.interactions$coefficients["age.cat70-79"],
                              outmod.interactions$coefficients["age.cat80-89"], outmod.interactions$coefficients["age.cat90-"])), matrix(0,nrow = 4,ncol=1),
                     matrix(c(outmod.interactions$coefficients["sex1:lowered.consc"], outmod.interactions$coefficients["age.cat70-79:lowered.consc"],
                              outmod.interactions$coefficients["age.cat80-89:lowered.consc"], outmod.interactions$coefficients["age.cat90-:lowered.consc"])), matrix(0,nrow = 4,ncol=1))

  data.new <-  RSdata
  data.new$af <- as.numeric(data.new$AF) - 1
  data.new$my <- data.new$af*data.new$lowered.consc
  medmod.order <- glm(lowered.consc ~ AF * (sex + age.cat), data = RSdata, family = binomial(link = "probit"))
  outmod.order <- glm(cf.3mo ~ AF + lowered.consc + my + lowered.consc*(sex + age.cat), data = data.new, family = binomial(link = "probit"))

  glm.sigma.thetabeta <- matrix(0, nrow = (2 + 2*4 + 4 + 4*4), ncol = (2 + 2*4 + 4 + 4*4))
  glm.sigma.thetabeta[21:30, 21:30] <- vcov(medmod.order)
  glm.sigma.thetabeta[c(1:8,13:16), c(1:8,13:16)] <- vcov(outmod.order)
  glm.sigma.thetabeta <- list(glm.sigma.thetabeta)

  expect_equivalent(sensmed.interactions$betas, glm.betas)
  expect_equivalent(sensmed.interactions$thetas, glm.thetas)
  expect_equivalent(sensmed.interactions$sigma.thetabeta, glm.sigma.thetabeta)

})

test_that("coefficients/covariances correctly stored for mediator and outcome models with interactions involving the exposure and/or mediator, one covariate", {

  medmod.int1cov <- glm(lowered.consc ~ AF * sex, data = RSdata, family = binomial(link = "probit"))
  outmod.int1cov <- glm(cf.3mo ~ sex * AF * lowered.consc, data = RSdata, family = binomial(link = "probit"))

  sensmed.int1cov <- sensmediation(med.model = medmod.int1cov, out.model = outmod.int1cov, exp.name = "AF1",
                                   med.name = "lowered.consc")
  sensmed.int1cov$betas
  sensmed.int1cov$thetas
  sensmed.int1cov$sigma.thetabeta

  glm.betas.int1cov <- list(medmod.int1cov$coefficients["(Intercept)"],medmod.int1cov$coefficients["AF1"],
                            matrix(medmod.int1cov$coefficients["sex1"]), matrix(medmod.int1cov$coefficients["AF1:sex1"]))
  glm.thetas.int1cov <- list(outmod.int1cov$coefficients["(Intercept)"],outmod.int1cov$coefficients["AF1"],
                             outmod.int1cov$coefficients["lowered.consc"], outmod.int1cov$coefficients["AF1:lowered.consc"],
                             matrix(outmod.int1cov$coefficients["sex1"]), matrix(outmod.int1cov$coefficients["sex1:AF1"]),
                             matrix(outmod.int1cov$coefficients["sex1:lowered.consc"]), matrix(outmod.int1cov$coefficients["sex1:AF1:lowered.consc"]))

  data.new <-  RSdata
  data.new$af <- as.numeric(data.new$AF) - 1
  data.new$my <- data.new$af*data.new$lowered.consc

  outmod.order.int1cov <- glm(cf.3mo ~ AF + lowered.consc + my + AF*sex + lowered.consc*sex + my*sex, data = data.new, family = binomial(link = "probit"))

  glm.sigma.thetabeta.int1cov <- matrix(0, nrow = (2 + 2*1 + 4 + 4*1), ncol = (2 + 2*1 + 4 + 4*1))
  glm.sigma.thetabeta.int1cov[9:12, 9:12] <- vcov(medmod.int1cov)
  glm.sigma.thetabeta.int1cov[c(1:8), c(1:8)] <- vcov(outmod.order.int1cov)
  glm.sigma.thetabeta.int1cov <- list(glm.sigma.thetabeta.int1cov)

  expect_equivalent(sensmed.int1cov$betas, glm.betas.int1cov)
  expect_equivalent(sensmed.int1cov$thetas, glm.thetas.int1cov)
  expect_equivalent(sensmed.int1cov$sigma.thetabeta, glm.sigma.thetabeta.int1cov)

})

test_that("effects correct for mediator and outcome models with interactions involving the exposure and/or mediator for Rho = 0 and Rho != 0", {

  medmod.interactions <- glm(lowered.consc ~ AF + sex + age.cat + AF*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))
  outmod.interactions <- glm(cf.3mo ~ sex + age.cat + AF * lowered.consc + lowered.consc*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))

  data.new <-  RSdata
  data.new$af <- as.numeric(data.new$AF) - 1
  data.new$my <- data.new$af*data.new$lowered.consc

  medmod.order <- glm(lowered.consc ~ AF * (sex + age.cat), data = RSdata, family = binomial(link = "probit"))
  outmod.order2 <- glm(cf.3mo ~ AF * lowered.consc + lowered.consc*(sex + age.cat), data = data.new, family = binomial(link = "probit"))

  sensmed.interactions0.1 <- sensmediation(med.model = medmod.interactions, out.model = outmod.interactions, exp.name = "AF1",
                                           med.name = "lowered.consc", Rho = c(0, 0.1))
  sensmed.order0.1 <- sensmediation(med.model = medmod.order, out.model = outmod.order2, exp.name = "AF1",
                                    med.name = "lowered.consc", Rho = c(0, 0.1))
  expect_equivalent(sensmed.order0.1$betas, sensmed.interactions0.1$betas)
  expect_equivalent(sensmed.order0.1$thetas, sensmed.interactions0.1$thetas)
  expect_equivalent(sensmed.order0.1$sigma.thetabeta, sensmed.interactions0.1$sigma.thetabeta)
  expect_equal(sensmed.order0.1$NIE, sensmed.interactions0.1$NIE)
  expect_equal(sensmed.order0.1$CI, sensmed.interactions0.1$CI)

})

test_that("errors are thrown for incorrectly specified exposure or mediator names (exp.name, med.name)", {

  medmod.1covariate <- glm(lowered.consc ~ AF + sex, data = RSdata, family = binomial(link = "probit"))
  outmod.1covariate <- glm(cf.3mo ~ AF + lowered.consc + sex, data = RSdata, family = binomial(link = "probit"))

  expect_error(sensmediation(med.model = medmod.1covariate, out.model = outmod.1covariate, exp.name = "AF",
                             med.name = "lowered.consc"), "The exposure is either missing")
  expect_error(sensmediation(med.model = medmod.1covariate, out.model = outmod.1covariate, exp.name = "AF1",
                             med.name = "lowered.consc1"), "The mediator is either missing")


})

test_that("alt.decomposition = TRUE is equal to the negative of switching the control and exposure values", {

  medmod.interactions <- glm(lowered.consc ~ AF + sex + age.cat + AF*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))
  outmod.interactions <- glm(cf.3mo ~ sex + age.cat + AF * lowered.consc + lowered.consc*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))

  sensmed.interactions <- sensmediation(med.model = medmod.interactions, out.model = outmod.interactions, exp.name = "AF1",
                                        med.name = "lowered.consc")

  sensmed.interactions.alt <- more.effects(sensmed.interactions, alt.decomposition = TRUE)
  sensmed.interactions.rev <- more.effects(sensmed.interactions, control.value = 1, exp.value = 0)
  expect_equivalent(sensmed.interactions.alt$NIE, -sensmed.interactions.rev$NIE)
  expect_equivalent(sensmed.interactions.alt$NDE, -sensmed.interactions.rev$NDE)


})

test_that("tests for conditional effects, effect estimates correct, error thrown for invalid covariate values, message for unused covariates", {

  medmod.interactions <- glm(lowered.consc ~ AF + sex + age.cat + AF*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))
  outmod.interactions <- glm(cf.3mo ~ sex + age.cat + AF * lowered.consc + lowered.consc*(age.cat + sex), data = RSdata, family = binomial(link = "probit"))

  sensmed.interactions <- sensmediation(med.model = medmod.interactions, out.model = outmod.interactions, exp.name = "AF1",
                                        med.name = "lowered.consc")

  sensmed.cond <- more.effects(sensmed.interactions, covariates = list(sex = 1, age.cat = "70-79"))

  b0 <- medmod.interactions$coefficients["(Intercept)"]
  b1 <- medmod.interactions$coefficients["AF1"]
  b2 <- matrix(c(medmod.interactions$coefficients["sex1"], medmod.interactions$coefficients["age.cat70-79"],
                 medmod.interactions$coefficients["age.cat80-89"], medmod.interactions$coefficients["age.cat90-"]))
  b3 <- matrix(c(medmod.interactions$coefficients["AF1:sex1"], medmod.interactions$coefficients["AF1:age.cat70-79"],
                 medmod.interactions$coefficients["AF1:age.cat80-89"], medmod.interactions$coefficients["AF1:age.cat90-"]))

  th0 <- outmod.interactions$coefficients["(Intercept)"]
  th1 <- outmod.interactions$coefficients["AF1"]
  th2 <- outmod.interactions$coefficients["lowered.consc"]
  th3 <- outmod.interactions$coefficients["AF1:lowered.consc"]
  th4 <- matrix(c(outmod.interactions$coefficients["sex1"], outmod.interactions$coefficients["age.cat70-79"],
                  outmod.interactions$coefficients["age.cat80-89"], outmod.interactions$coefficients["age.cat90-"]))
  th6 <- matrix(c(outmod.interactions$coefficients["sex1:lowered.consc"], outmod.interactions$coefficients["age.cat70-79:lowered.consc"],
                  outmod.interactions$coefficients["age.cat80-89:lowered.consc"], outmod.interactions$coefficients["age.cat90-:lowered.consc"]))
  x <- matrix(c(1,1,0,0), nrow = 1)


  nie <- (stats::pnorm(th0 + th2 + th1 + th3 + x%*%(th4 + th6))- stats::pnorm(th0 + th1 + x%*%th4))*(stats::pnorm(b0 + b1 + x%*%(b2 + b3))- stats::pnorm(b0 + x%*%b2))
  nde <- (stats::pnorm(th0 + th1 + x%*%th4) - stats::pnorm(th0 + x%*%th4))*(1 - stats::pnorm(b0 + x%*%b2 ) )  +
    (stats::pnorm(th0 + th2 + th1 + th3 + x%*%(th4 + th6)) - stats::pnorm(th0 + th2 + x%*%(th4 + th6)))*stats::pnorm(b0 + x%*%b2 )

  expect_equivalent(sensmed.cond$NIE, nie)
  expect_equivalent(sensmed.cond$NDE, nde)

  sensmed.cond2 <- more.effects(sensmed.interactions, covariates = list(age.cat = "70-79"))

  x2 <- matrix(c(as.numeric(RSdata$sex)-1, rep(1,1000), rep(0,1000), rep(0,1000)), nrow = 1000)

  nie2 <- matrix(mean((stats::pnorm(th0 + th2 + th1 + th3 + x2%*%(th4 + th6))- stats::pnorm(th0 + th1 + x2%*%th4))*(stats::pnorm(b0 + b1 + x2%*%(b2 + b3))- stats::pnorm(b0 + x2%*%b2))))
  nde2 <- matrix(mean((stats::pnorm(th0 + th1 + x2%*%th4) - stats::pnorm(th0 + x2%*%th4))*(1 - stats::pnorm(b0 + x2%*%b2 ) )  +
                        (stats::pnorm(th0 + th2 + th1 + th3 + x2%*%(th4 + th6)) - stats::pnorm(th0 + th2 + x2%*%(th4 + th6)))*stats::pnorm(b0 + x2%*%b2 ) ))

  expect_equivalent(sensmed.cond2$NIE, nie2)
  expect_equivalent(sensmed.cond2$NDE, nde2)

  expect_error(more.effects(sensmed.interactions, covariates = list(sex = 2, age.cat = "70-79")), "is not a level of" )

  expect_message(more.effects(sensmed.interactions, covariates = list(sex = 1, age.cat = "70-79", bla = 2, apa = 3)), "Note: bla, apa not found")

  medmod.int1cov <- glm(lowered.consc ~ AF * sex, data = RSdata, family = binomial(link = "probit"))
  outmod.int1cov <- glm(cf.3mo ~ sex * AF * lowered.consc, data = RSdata, family = binomial(link = "probit"))

  sensmed.int1cov <- sensmediation(med.model = medmod.int1cov, out.model = outmod.int1cov, exp.name = "AF1",
                                   med.name = "lowered.consc")

  sensmed.cond3 <- more.effects(sensmed.int1cov, covariates = list(sex = 1))
  expect_equivalent(round(sensmed.cond3$NIE,8), matrix(0.01596528))

  data.new <-  RSdata
  data.new$age.char <- as.character(data.new$age.cat)
  medmod.agechar <- glm(lowered.consc ~ AF + sex + age.char, data = data.new, family = binomial(link = "probit"))
  outmod.agechar <- glm(cf.3mo ~ AF + lowered.consc + sex + age.char, data = data.new, family = binomial(link = "probit"))

  sensmed.agechar <- sensmediation(med.model = medmod.agechar, out.model = outmod.agechar, exp.name = "AF1",
                                   med.name = "lowered.consc")

  expect_error(more.effects(sensmed.agechar, covariates = list(sex = 1, age.char = "70-79")), "may not be of class" )

})




