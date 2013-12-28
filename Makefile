export GOBIN=./bin

fmt:
	gofmt -w cmd/*/*.go
	colcheck cmd/*/*.go

push:
	git push origin master
	git push tufts master
	git push github master

