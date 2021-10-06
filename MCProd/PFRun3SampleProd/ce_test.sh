#!/bin/bash

export X509_USER_PROXY=$PWD/x509up

date
hostname
env
voms-proxy-info

#sleep 30
