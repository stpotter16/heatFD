#!/usr/bin/env bats

executable="../src/heatFD"

@test "verify heatFD executable exists" {
	run ls $executable
	[ "$status" -eq 0 ]
}
