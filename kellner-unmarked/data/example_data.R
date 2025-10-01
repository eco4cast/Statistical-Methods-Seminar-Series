library(unmarked)
data(crossbill)

# Detection/nondetection
y99 <- crossbill[,c("det991", "det992", "det993")]
y00 <- crossbill[,c("det001", "det002", "det003")]
y01 <- crossbill[,c("det011", "det012", "det013")]
colnames(y99) <- colnames(y00) <- colnames(y01) <- c("occ1", "occ2", "occ3")
y <- rbind(y99, y00, y01)

site_covs <- crossbill[,c("id", "forest", "ele")]
# Replicate 3 times
site_covs <- rbind(site_covs, site_covs, site_covs)

date99 <- crossbill[,c("date991", "date992", "date993")]
date00 <- crossbill[,c("date001", "date002", "date003")]
date01 <- crossbill[,c("date011", "date012", "date013")]
colnames(date99) <- colnames(date00) <- colnames(date01) <- c("date1", "date2", "date3")
dates <- rbind(date99, date00, date01)

write.csv(y, "y.csv", row.names=FALSE)
write.csv(site_covs, "site_covs.csv", row.names=FALSE)
write.csv(dates, "dates.csv", row.names = FALSE)
