maizeSNPDB
========

The genotypes of 1210 maize lines at 35,370,939 SNP sites are stored using Sparse Matrices in R.  

*****

#	Use maizeSNPDB online

maizeSNPDB is deployed at <a href="http://150.109.59.144:3838/MaizeSNPDB/" target="_blank">http://150.109.59.144:3838/MaizeSNPDB/</a> for online use.  
maizeSNPDB is idle until you activate it by accessing the URL. So it may take some time when you access this URL for the first time. Once it was activated, maizeSNPDB could be used smoothly and easily.  

*****

#	<font color="red">Launch maizeSNPDB directly from R and GitHub (preferred approach)</font>

User can choose to run maizeSNPDB installed locally for a more preferable experience.

**Step 1: Install R and RStudio**

Before running the app you will need to have R and RStudio installed (tested with R 3.4.4 and RStudio 1.1.442).  
Please check CRAN (<a href="https://cran.r-project.org/" target="_blank">https://cran.r-project.org/</a>) for the installation of R.  
Please check <a href="https://www.rstudio.com/" target="_blank">https://www.rstudio.com/</a> for the installation of RStudio.  

**Step 2: Install the R Shiny package and other packages required by maizeSNPDB**

Start an R session using RStudio and run these lines:  
```
# try an http CRAN mirror if https CRAN mirror doesn't work  
install.packages("shiny")  
install.packages("shinyBS")  
install.packages("shinythemes")  
install.packages("shinycssloaders")  
install.packages("plotly")  
install.packages("foreach")  
install.packages("ape")  
install.packages("pegas")  
install.packages("dplyr")  
install.packages("ggmap")  
install.packages("tidyr")  
install.packages("gridExtra")  
install.packages("genetics")  
# try http:// if https:// URLs are not supported   
source("https://bioconductor.org/biocLite.R")  
biocLite("IRanges")
biocLite("snpStats")
biocLite("chopsticks")  
biocLite("ggtree")  
# try an http CRAN mirror if https CRAN mirror doesn't work  
install.packages("LDheatmap")
# install shinysky  
if (require(devtools)) install.packages("devtools")  
devtools::install_github("venyao/ShinySky")  
```

**Step 3: Start the app**  

Start an R session using RStudio and run these lines:  
```
library(shiny)  
runGitHub("maizeSNPDB", "venyao", launch.browser = TRUE)  
```
This command would take some time as it will download the maizeSNPDB database from GitHub to the disk of your local computer (check the directory path using the function `getwd()` in R).   

<br>
Alternatively, you can download the maizeSNPDB database from Jianguoyun (https://www.jianguoyun.com/p/DXP-bAQQzqnhBRiu8Vg) or GitHub (https://github.com/venyao/maizeSNPDB) to a directory (for example "E:/apps/") of your local computer using the web browser or other tools.   

<br>
<img src="maizeSNPDB.png" width="890"/>  
<br>

Then start an R session using RStudio and run these lines:  
```
library(shiny)  
runApp("E:/apps/maizeSNPDB", launch.browser = TRUE)  
# The first parameter of runApp should be the directory that contains the scripts server.R and ui.R of maizeSNPDB.  
```

Your web browser will open the app.

*****

#	Deploy maizeSNPDB on local or web Linux server

**Step 1: Install R**  

Please check CRAN (<a href="https://cran.r-project.org/" target="_blank">https://cran.r-project.org/</a>) for the installation of R.

**Step 2: Install the R Shiny package and other packages required by maizeSNPDB**  

Start an R session and run these lines in R:  
```
# try an http CRAN mirror if https CRAN mirror doesn't work  
install.packages("shiny")  
install.packages("shinyBS")  
install.packages("shinythemes")  
install.packages("shinycssloaders")  
install.packages("plotly")  
install.packages("foreach")  
install.packages("ape")  
install.packages("pegas")  
install.packages("dplyr")  
install.packages("ggmap")  
install.packages("tidyr")  
install.packages("gridExtra")  
install.packages("genetics")  
# try http:// if https:// URLs are not supported   
source("https://bioconductor.org/biocLite.R")  
biocLite("IRanges")
biocLite("snpStats")
biocLite("chopsticks")  
biocLite("ggtree")  
# try an http CRAN mirror if https CRAN mirror doesn't work  
install.packages("LDheatmap")  
# install shinysky  
if (require(devtools)) install.packages("devtools")  
devtools::install_github("venyao/ShinySky")  
```

For more information, please check the following pages:  
<a href="https://cran.r-project.org/web/packages/shiny/index.html" target="_blank">https://cran.r-project.org/web/packages/shiny/index.html</a>  
<a href="https://github.com/rstudio/shiny" target="_blank">https://github.com/rstudio/shiny</a>  
<a href="https://shiny.rstudio.com/" target="_blank">https://shiny.rstudio.com/</a>  

**Step 3: Install Shiny-Server**

Please check the following pages for the installation of shiny-server.  
<a href="https://www.rstudio.com/products/shiny/download-server/" target="_blank">https://www.rstudio.com/products/shiny/download-server/</a>  
<a href="https://github.com/rstudio/shiny-server/wiki/Building-Shiny-Server-from-Source" target="_blank">https://github.com/rstudio/shiny-server/wiki/Building-Shiny-Server-from-Source</a>  

**Step 4: Upload files of maizeSNPDB**

Put the directory containing the code and data of maizeSNPDB to /srv/shiny-server.  

**Step 5: Configure shiny server (/etc/shiny-server/shiny-server.conf)**

```
# Define the user to spawn R Shiny processes
run_as shiny;

# Define a top-level server which will listen on a port
server {  
  # Use port 3838  
  listen 3838;  
  # Define the location available at the base URL  
  location /maizesnp {  
    # Directory containing the code and data of maizeSNPDB  
    app_dir /srv/shiny-server/maizeSNPDB;  
    # Directory to store the log files  
    log_dir /var/log/shiny-server;  
  }  
}  
```

**Step 6: Change the owner of the maizeSNPDB directory**

```
$ chown -R shiny /srv/shiny-server/maizeSNPDB  
```

**Step 7: Start Shiny-Server**

```
$ start shiny-server  
```

Now, the maizeSNPDB app is available at http://IPAddressOfTheServer:3838/maizeSNPDB/.  

