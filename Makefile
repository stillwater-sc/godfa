

api-server:
	go install github.com/stillwater-sc/godfa/cmd/dfa-api-server


indexspace:
	go test ./...
