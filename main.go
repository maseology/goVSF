package main

import (
	"fmt"
	"math"

	. "github.com/maseology/goHydro/profile"
	. "github.com/maseology/goHydro/richards1d"
)

func main() {
	simLength := 5.0 // hours
	maxTimeStep := 3600.0
	initSe := 0.3

	var p Profile
	p.New()
	// fmt.Printf("soil profile {ts, tr, ks, he, b} = %v\n", *p.P[1])

	run(p, simLength, maxTimeStep, initSe)
}

func run(p Profile, simLength, maxTimeStep, initSe float64) {
	// control
	endTime := simLength * 3600.0 // sec
	dt := maxTimeStep / 10.0
	time := 0.0
	sumInfiltration := 0.0
	totiter := 0
	// boundary conditions
	ubPotential := p.P[1].He
	isFreeDrainage := true
	// state variables
	var t ProfileState
	t.InitializeWater(p, initSe)

	for time < endTime {
		dt = math.Min(dt, endTime-time)
		ok, inneriter, infil := t.CellCentFiniteVolWater(p, dt, ubPotential, isFreeDrainage)
		totiter += inneriter
		if ok {
			for i := 0; i <= NsubLay+1; i++ {
				t.Tl[i] = t.T[i]
			}
			sumInfiltration += infil * dt
			time += dt
			fmt.Printf(" time = %d\tdt = %.2f\titer = %d\tinfil = %.3f\n", int(time), dt, inneriter, sumInfiltration)

			if float64(inneriter)/float64(MaxIter) < 0.1 {
				dt = math.Min(dt*2.0, maxTimeStep)
			}
		} else {
			fmt.Printf("dt = %.3f\tNo convergence\n", dt)
			dt = math.Max(dt/2.0, 1.0)
			for i := 0; i <= NsubLay+1; i++ {
				t.T[i] = t.Tl[i]
				t.Psi[i] = t.PM[i].GetPsi(t.T[i])
			}
		}
	}
	fmt.Printf("number of iterations per hour: %.1f\n\n", float64(totiter)/simLength)
}
