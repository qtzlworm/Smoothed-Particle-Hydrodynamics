/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
/*
 * owPhycisTest.cpp
 *
 *  Created on: Apr 29, 2014
 *      Author: Sergey Khayrulin
 *      email: s.khayrulin@gmail.com
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>        // std::abs

#include "owPhysicTest.h"
#include "../owPhysicsFluidSimulator.h"

float calcPotentialEnergy( owConfigProrerty *, float * );
float calcKineticEnergy( owConfigProrerty *, float *, float * );
float get_len( float * );
void get_position(owPhysicsFluidSimulator * , float * , float * , int );
float get_dist(float *, float *);
float gravity = 9.81f;

/*******************************************************
 * CONSERVATION ENERGY TEST DESCRIPTION
 * TOTAL ENERGY OF SYSTEM SHOULD BE CONSTANT ALL TIME
 * E = mv^2/2 + mgh
 * *****************************************************/
void test_energy_conservation(){
	owHelper::path = "./configuration/test/";
	owHelper::suffix = "_liquid";
	owHelper * helper = new owHelper();
	owPhysicsFluidSimulator * fluid_simulation = new owPhysicsFluidSimulator(helper);
	float total_energy = 0.f;
	float kinetic_energy = 0.f;
	float potential_energy = 0.f;
	float * p_buffer;
	float * v_buffer;
	std::vector<float> energy_evolution_total;
	std::vector<float> energy_evolution_potential;
	std::vector<float> energy_evolution_kinetic;
	int counter = 0;
	std::cout << "===================" << "CONSERVATION ENERGY TEST START" << "========================" << std::endl;
	while(1){
		p_buffer = fluid_simulation->getPosition_cpp();
		v_buffer = fluid_simulation->getvelocity_cpp();
		potential_energy = calcPotentialEnergy(fluid_simulation->getConfig(),p_buffer);
		kinetic_energy = calcKineticEnergy(fluid_simulation->getConfig(),v_buffer,p_buffer);
		total_energy = kinetic_energy + potential_energy;
		energy_evolution_total.push_back(total_energy);
		energy_evolution_kinetic.push_back(kinetic_energy);
		energy_evolution_potential.push_back(potential_energy);
		fluid_simulation->simulationStep();
		if(counter == 10000)
			break;
		counter++;
	}
	owHelper::log_buffer(&energy_evolution_total[0], 1, energy_evolution_total.size(), "./logs/total_energy_distrib.txt");
	owHelper::log_buffer(&energy_evolution_kinetic[0], 1, energy_evolution_kinetic.size(), "./logs/kinetic_energy_distrib.txt");
	owHelper::log_buffer(&energy_evolution_potential[0], 1, energy_evolution_potential.size(), "./logs/potential_energy_distrib.txt");
	std::cout << "===================" << "CONSERVATION ENERGY TEST END  " << "========================" << std::endl;
	delete fluid_simulation;
}

void test_gravity(){
	owHelper::path = "./configuration/test/";
	owHelper::suffix = "_elastic_one";
	owHelper * helper = new owHelper();
	owPhysicsFluidSimulator * fluid_simulation = new owPhysicsFluidSimulator(helper);
	float * initial_position = new float[4];
	float * initial_velocity = new float[4];
	float * current_position = new float[4];
	float * current_velocity = new float[4];
	float * p_buffer;
	float * v_buffer;
	const int id = 0;
	int counter = 0;
	const int totalNumberOfIteration = 1000;
	float * result = new float[ totalNumberOfIteration * 2 ];
	get_position(fluid_simulation,initial_position,initial_velocity,id);
	result[0] = timeStep * counter;
	result[1] = get_dist(initial_position,initial_position) * simulationScale;
	counter++;
	std::cout << "===================" << "CONSERVATION GRAVITY TEST START" << "========================" << std::endl;
	while(1){
		fluid_simulation->simulationStep();
		get_position(fluid_simulation,current_position,current_velocity,id);
		if(counter == totalNumberOfIteration)
			break;
		result[counter * 2 + 0] = timeStep * counter;
		result[counter * 2 + 1] = get_dist(initial_position,current_position) * simulationScale;
		counter++;
	}
	owHelper::log_buffer(result, 2, totalNumberOfIteration, "./logs/gravity_test_distrib.txt");
	std::cout << "===================" << "CONSERVATION GRAVITY TEST END  " << "========================" << std::endl;
	delete fluid_simulation;
}
float calcPotentialEnergy(owConfigProrerty * config, float * p_buffer){
	float e = 0.f;
	float l = 0.f;
	for(int i=0;i<config->getParticleCount();i++){
		if((int)(p_buffer[4 * i + 3]) != BOUNDARY_PARTICLE){
			l = (p_buffer[4 * i + 1] <= r0) ? 0.f : p_buffer[4 * i + 1] * simulationScale;
			e += std::abs(l) * gravity * mass; //Y - coordinate is a h
		}
	}
	return e;
}
float calcKineticEnergy(owConfigProrerty * config, float * v_buffer, float * p_buffer){
	float e = 0.f;
	for(int i=0;i < config->getParticleCount();i++){
		if((int)(p_buffer[4 * i + 3]) != BOUNDARY_PARTICLE){
			e += mass * pow(get_len(v_buffer + 4 * i + 0),2.0f)/2.0f; //
		}
	}
	return e;
}
float get_dist(float * from, float * to){
	return sqrt ( pow(to[0] - from[0], 2.f) + pow(to[1] - from[1], 2.f) + pow(to[2] - from[2], 2.f) );
}
float get_len(float * v){
	return sqrt( pow(v[0],2.0f) + pow(v[1],2.0f) + pow(v[2],2.0f));
}
void get_position(owPhysicsFluidSimulator * fluid_simulation, float * init_p, float * init_v, int id){
	float * p_b;
	float * v_b;
	p_b = fluid_simulation->getPosition_cpp();
	v_b = fluid_simulation->getvelocity_cpp();
	init_p[0] = p_b[ id * 4 + 0 ];
	init_p[1] = p_b[ id * 4 + 1 ];
	init_p[2] = p_b[ id * 4 + 2 ];
	init_p[3] = p_b[ id * 4 + 3 ];

	init_v[0] = v_b[ id * 4 + 0 ];
	init_v[1] = v_b[ id * 4 + 1 ];
	init_v[2] = v_b[ id * 4 + 2 ];
	init_v[3] = v_b[ id * 4 + 3 ];
	//delete v_b;
	//delete p_b;
}
