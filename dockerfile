FROM rstan/covid

# Install 
RUN install2.r --error \
        pracma 
    

COPY /root/.shh/ /root/.shh/


# docker run -d -i -t -p 8787:8787 -u root -v /covidmodeldata:/covidmodeldata -e DISABLE_AUTH=true -e ROOT=TRUE -e USER=cades rstan/covid-master


# docker run --rm -it --volume /home/cades:/home/rstudio/covid-model/cades --env USER=cades --env PASSWORD=password --env ROOT=TRUE rstan/covid /bin/bash