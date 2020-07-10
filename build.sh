#!/usr/bin/env bash

env GOOS=darwin GOARCH=amd64 go build -o REDItools2_darwin_amd64 .
env GOOS=linux GOARCH=amd64 go build -o REDItools2_linux_amd64 .
env GOOS=windows GOARCH=amd64 go build -o REDItools2_win_amd64.exe .

