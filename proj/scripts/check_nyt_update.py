#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check the nytimes/covid-19-data git repository for updated file."""
import time
import subprocess
from github import Github

# define the command to run
CMD = ['Rscript', '/home/cades/covid-model/run_model.R']

url = 'https://github.com/nytimes/covid-19-data/'
repo = 'nytimes/covid-19-data'
county_csv = 'us-counties.csv'


def retrieve_date():
    """Retrieve the date from the nytimes Github."""
    g = Github()
    county_date = g.get_repo('nytimes/covid-19-data'
                             ).get_commits(
                                 path=county_csv)[0].commit.committer.date
    return county_date


if __name__ == "__main__":
    county_date = retrieve_date()
    while True:
        if county_date != retrieve_date():
            process = subprocess.Popen(CMD, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            county_date = retrieve_date()
        time.sleep(3600)
