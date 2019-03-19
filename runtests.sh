#!/bin/bash

py.test --pyargs speedcom --cov-report term-missing --cov=speedcom --cov-config .coveragerc 
