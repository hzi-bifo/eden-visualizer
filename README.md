# eden-visualizer

online-hosted: https://pmuench.shinyapps.io/eden-visualizer/

when using eden, this application will run inside a docker container and can be accessed using localhost

# docker version

you can run it inside a Docker container by adding this to the Dockerfile

```
RUN wget https://rawgit.com/philippmuench/eden-visualizer/master/bundle.tar.gz -O /srv/shiny-server/bundle.tar.gz &&\
  tar -xvzf /srv/shiny-server/bundle.tar.gz --directory=/srv/shiny-server/
```
