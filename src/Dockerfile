FROM ubuntu
# Make sure local repositories are up-to-date, for apt to find needed packages
RUN apt-get update
# -y option avoids the terminal prompting confirmation, which we wont be able to answer to as the build is not interactive
RUN apt-get install python3 python3-pip -y
RUN python3 -m pip install --extra-index-url https://test.pypi.org/simple/ HandsOnDeployAPackage
WORKDIR /root/home

