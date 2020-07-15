#!/usr/bin/env bash

git checkout master
env GOOS=darwin GOARCH=amd64 go build -o REDItools2_l_darwin_amd64 .
env GOOS=linux GOARCH=amd64 go build -o REDItools2_l_linux_amd64 .
env GOOS=windows GOARCH=amd64 go build -o REDItools2_l_win_amd64.exe .


git checkout dev
env GOOS=darwin GOARCH=amd64 go build -o REDItools2_h_darwin_amd64 .
env GOOS=linux GOARCH=amd64 go build -o REDItools2_h_linux_amd64 .
env GOOS=windows GOARCH=amd64 go build -o REDItools2_h_win_amd64.exe .

git checkout master
