
const fs = require('fs');

/*
  I've wrapped Makoto Matsumoto and Takuji Nishimura's code in a namespace
  so it's better encapsulated. Now you can have multiple random number generators
  and they won't stomp all over eachother's state.
  
  If you want to use this as a substitute for Math.random(), use the random()
  method like so:
  
  var m = new MersenneTwister();
  var randomNumber = m.random();
  
  You can also call the other genrand_{foo}() methods on the instance.

  If you want to use a specific seed in order to get a repeatable random
  sequence, pass an integer into the constructor:

  var m = new MersenneTwister(123);

  and that will always produce the same random sequence.

  Sean McCullough (banksean@gmail.com)
*/

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
 
   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).
 
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          
 
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
 
     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
 
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
 
     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.
 
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 
   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

var MersenneTwister = function(seed) {
  if (seed == undefined) {
    seed = new Date().getTime();
  } 
  /* Period parameters */  
  this.N = 624;
  this.M = 397;
  this.MATRIX_A = 0x9908b0df;   /* constant vector a */
  this.UPPER_MASK = 0x80000000; /* most significant w-r bits */
  this.LOWER_MASK = 0x7fffffff; /* least significant r bits */
 
  this.mt = new Array(this.N); /* the array for the state vector */
  this.mti=this.N+1; /* mti==N+1 means mt[N] is not initialized */

  this.init_genrand(seed);
}  
 
/* initializes mt[N] with a seed */
MersenneTwister.prototype.init_genrand = function(s) {
  this.mt[0] = s >>> 0;
  for (this.mti=1; this.mti<this.N; this.mti++) {
      var s = this.mt[this.mti-1] ^ (this.mt[this.mti-1] >>> 30);
   this.mt[this.mti] = (((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253)
  + this.mti;
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      this.mt[this.mti] >>>= 0;
      /* for >32 bit machines */
  }
}
 
/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
MersenneTwister.prototype.init_by_array = function(init_key, key_length) {
  var i, j, k;
  this.init_genrand(19650218);
  i=1; j=0;
  k = (this.N>key_length ? this.N : key_length);
  for (; k; k--) {
    var s = this.mt[i-1] ^ (this.mt[i-1] >>> 30)
    this.mt[i] = (this.mt[i] ^ (((((s & 0xffff0000) >>> 16) * 1664525) << 16) + ((s & 0x0000ffff) * 1664525)))
      + init_key[j] + j; /* non linear */
    this.mt[i] >>>= 0; /* for WORDSIZE > 32 machines */
    i++; j++;
    if (i>=this.N) { this.mt[0] = this.mt[this.N-1]; i=1; }
    if (j>=key_length) j=0;
  }
  for (k=this.N-1; k; k--) {
    var s = this.mt[i-1] ^ (this.mt[i-1] >>> 30);
    this.mt[i] = (this.mt[i] ^ (((((s & 0xffff0000) >>> 16) * 1566083941) << 16) + (s & 0x0000ffff) * 1566083941))
      - i; /* non linear */
    this.mt[i] >>>= 0; /* for WORDSIZE > 32 machines */
    i++;
    if (i>=this.N) { this.mt[0] = this.mt[this.N-1]; i=1; }
  }

  this.mt[0] = 0x80000000; /* MSB is 1; assuring non-zero initial array */ 
}
 
/* generates a random number on [0,0xffffffff]-interval */
MersenneTwister.prototype.genrand_int32 = function() {
  var y;
  var mag01 = new Array(0x0, this.MATRIX_A);
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (this.mti >= this.N) { /* generate N words at one time */
    var kk;

    if (this.mti == this.N+1)   /* if init_genrand() has not been called, */
      this.init_genrand(5489); /* a default initial seed is used */

    for (kk=0;kk<this.N-this.M;kk++) {
      y = (this.mt[kk]&this.UPPER_MASK)|(this.mt[kk+1]&this.LOWER_MASK);
      this.mt[kk] = this.mt[kk+this.M] ^ (y >>> 1) ^ mag01[y & 0x1];
    }
    for (;kk<this.N-1;kk++) {
      y = (this.mt[kk]&this.UPPER_MASK)|(this.mt[kk+1]&this.LOWER_MASK);
      this.mt[kk] = this.mt[kk+(this.M-this.N)] ^ (y >>> 1) ^ mag01[y & 0x1];
    }
    y = (this.mt[this.N-1]&this.UPPER_MASK)|(this.mt[0]&this.LOWER_MASK);
    this.mt[this.N-1] = this.mt[this.M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

    this.mti = 0;
  }

  y = this.mt[this.mti++];

  /* Tempering */
  y ^= (y >>> 11);
  y ^= (y << 7) & 0x9d2c5680;
  y ^= (y << 15) & 0xefc60000;
  y ^= (y >>> 18);

  return y >>> 0;
}
 
/* generates a random number on [0,0x7fffffff]-interval */
MersenneTwister.prototype.genrand_int31 = function() {
  return (this.genrand_int32()>>>1);
}
 
/* generates a random number on [0,1]-real-interval */
MersenneTwister.prototype.genrand_real1 = function() {
  return this.genrand_int32()*(1.0/4294967295.0); 
  /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
MersenneTwister.prototype.random = function() {
  return this.genrand_int32()*(1.0/4294967296.0); 
  /* divided by 2^32 */
}
 
/* generates a random number on (0,1)-real-interval */
MersenneTwister.prototype.genrand_real3 = function() {
  return (this.genrand_int32() + 0.5)*(1.0/4294967296.0); 
  /* divided by 2^32 */
}
 
/* generates a random number on [0,1) with 53-bit resolution*/
MersenneTwister.prototype.genrand_res53 = function() { 
  var a=this.genrand_int32()>>>5, b=this.genrand_int32()>>>6; 
  return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

/* These real versions are due to Isaku Wada, 2002/01/09 added */


function generateKenken (settings) {
	if (!settings.seed) settings.seed = new Date().getTime()
    
    var kenken = new Kenken(settings)
    // renderKenken(kenken)
    // location.hash = encodeOptions(settings)
    return kenken
}

// A class for the ken ken board
function Kenken (settings) {
    var size = this.size = settings.size
    this.settings = settings
    this.board = []
	this.minGroupSize = 1
	this.defaultMaxGroupSize = settings.maxGroupSize
	this.maxGroupSize = undefined
	this.cellGroups = []
	// Determine operation based on checkboxes from webpage
	this.operations = [new SingleCell()]
	if(settings.operations.addition) {
		this.operations.push(new Addition())
	}
	if(settings.operations.subtraction) {
		this.operations.push(new Subtraction())
	}
	if(settings.operations.multiplication) {
		this.operations.push(new Multiplication())
	}
	if(settings.operations.division) {
		this.operations.push(new Division())
	}
	if(settings.operations.max) {
		this.operations.push(new Maximum())
	}
	if(settings.operations.min) {
		this.operations.push(new Minimum())
	}
	if(settings.operations.range) {
		this.operations.push(new Range())
	}
	if(settings.operations.mod) {
		this.operations.push(new Mod())
	}
	if(settings.operations.avg) {
		this.operations.push(new Average())
	}
	if(settings.operations.par) {
		this.operations.push(new Parity())
	}
	if(settings.operations.gcd) {
		this.operations.push(new Gcd())
	}
	
	//Determine if the default maxGroupSize needs to be made smaller due to the operations selected
	for(var i = 0; i < this.operations.length; i++) {
		var maxOperationSize = this.operations[i].maxCells
		if(maxOperationSize == undefined) {
			// There is an operation with no maximum cell group size, so use the deafult max size (or the size from the webpage)
			this.maxGroupSize = this.defaultMaxGroupSize
			// Stop the loop
			i = this.operations.length
		} else if(this.maxGroupSize == undefined) {
			this.maxGroupSize = maxOperationSize
		} else if(maxOperationSize > this.maxGroupSize) {
			this.maxGroupSize = maxOperationSize
		}
	}
	// If the max group size determined on the webpage is smaller than the max from the operations, make this the max group size
	if(this.maxGroupSize > this.defaultMaxGroupSize) {
		this.maxGroupSize = this.defaultMaxGroupSize
	}
	
	this.seed = new MersenneTwister(settings.seed)
	
	var builderArray = shuffledNumberArray(size, this.seed)
    
    for (var x = 0; x < size; x++) {
        this.board[x] = []
        for (var y = 0; y < size; y++) {
			// Fills the board with cells, where each row of cells have the values
			// of the builderArray, cyclically shifted to the left by x (the row number)
            this.board[x][y] = new Cell(this, builderArray[(x+y)%size])
        }
    }
	// Shuffle the board
	shuffleBoard(size,this.board, this.seed)
	
	// Assign cells in the now shuffled board x and y values
	for(var x = 0; x < size; x++) {
		for(var y = 0; y < size; y++) {
			this.board[x][y].x = x;
			this.board[x][y].y = y;
		}
	}
	
	// Create the cell groups
	var groupID = 1
	for(var x = 0; x < size; x++) {
		for(var y = 0; y < size; y++) {
			if(this.board[x][y].cellGroup == undefined) {
				// Generate a random integer in the range [minGroupSize,maxGroupSize] for the size of the group
				var groupSize = Math.floor((this.maxGroupSize-this.minGroupSize+1)*this.seed.random()+(this.minGroupSize))
				// Create the new CellGroup object
				var newCellGroup = new CellGroup(this, this.board[x][y], groupID)
				// Grow the new cell group groupSize-1 times, only if the groupSize is not one (since it already has size one)
				if(groupSize != 1) {
					for(var m = 0; m < groupSize-1; m++) {
						newCellGroup.grow()
					}
				}
				
				// Shuffle the operations array to get a good spread of the operations, prevents "overshadowing" where having
				//  + before x picked + more often due to how the algorithm below selects the oepration
				shuffleInputArray(this.operations,this.seed)
				
				// Start at the beginning of this now random order operations array
				var randomOperationStart = 0 //Math.floor(this.operations.length*this.seed.random())
				var randomOperation = randomOperationStart
				// Runs until valid operation is found, should never infinite loop as there should always be a valid operation
				var foundOperation = false
				while(foundOperation == false) {
					if(this.operations[randomOperation].operation(newCellGroup.getAllValues()) != false) {
						// A valid operation was found, set the operation text for the group
						newCellGroup.operationDescription = this.operations[randomOperation].symbol + this.operations[randomOperation].operation(newCellGroup.getAllValues())
						foundOperation = true
					}
					else {
						// Not a valid operation, move on to the next options
						randomOperation = (randomOperation + 1) % this.operations.length
						if(randomOperation == randomOperationStart) {
							console.log('no valid operation was found')
							break
						}
					}
				}

				// Add the new cell group to the kenken cell group member variable
				this.cellGroups.push(newCellGroup)
				
				groupID = groupID + 1
			}
		}
	}
}

function shuffleBoard (size,board,seed) {
	// Swap two columns and then two rows. Do this 'size' times to get a decent mix up of the board.
	for (var i = 0; i < size; i++) {
		// Generate two random integers in the range [0,size)
		var column1 = Math.floor(size*seed.random())
		var column2 = Math.floor(size*seed.random())
		// Swap the two columns
		for(var j = 0; j < size; j++) {
			var tempCell = board[j][column1]
			board[j][column1] = board[j][column2]
			board[j][column2] = tempCell
		}
		
		// Generate two random integers in the range [0,size)
		var row1 = Math.floor(size*seed.random())
		var row2 = Math.floor(size*seed.random())
		// Swap the two rows
		for(var j = 0; j < size; j++) {
			var tempCell = board[row1][j]
			board[row1][j] = board[row2][j]
			board[row2][j] = tempCell
		}
	}
}

// A class for a single grouping of cells on the ken ken board
function CellGroup (kenken, cell, id) {
	// The ken ken board this group belongs to (is this necessary?)
	this.kenken = kenken
	// The id of this group of cells
	this.groupID = id
	// The array that will hold the cells in this group, starting with the initial cell
	this.cells = [cell]
	// Set the cellGroup of the initial cell to this CellGroup Object
	cell.setCellGroup(this)
	// The current size of the cell group
	this.currentSize = 1
	// The variable for the string showing the operation and result this cell should have
	this.operationDescription = undefined
}

// Return the values of all cells in the group as an array
CellGroup.prototype.getAllValues = function() {
	var returnArray = []
	for(var i = 0; i < this.cells.length; i++) {
		returnArray.push(this.cells[i].value)
	}
	return returnArray
}

// Grow the cell group up to maximum size, or smaller if board is not big enough
// Returns true if growing was successful, false if it was unsuccessful
CellGroup.prototype.grow = function() {
	// Generate a random integer in range [0,cells.length-1] for which cell we should attempt to grow at first
	var startingCellNumber = Math.floor(this.cells.length*this.kenken.seed.random())
	var cellNum = startingCellNumber
	while(true) {
		var cellToGrowFrom = this.cells[cellNum]
		// Get the array of neighbors of this cell
		var cellNeighbors = cellToGrowFrom.getNeighbors()
		//Generate a random integer in range [0,cellNeighbors.length-1]
		var neighborCellNum = Math.floor(cellNeighbors.length*this.kenken.seed.random())
		// Go through each neighbor. If one is valid, make it the next cell in this group.
		for(var i = 0; i < cellNeighbors.length; i++) {
			var neighborCell = cellNeighbors[((i+neighborCellNum)%cellNeighbors.length)]
			if(neighborCell.cellGroup == undefined) {
				this.cells.push(neighborCell)
				neighborCell.setCellGroup(this)
				this.currentSize = this.currentSize + 1 // Increase the current size count of the cell group
				return true
			}
		}
		
		// If all the neighbors were invalid, try the next cell in the list
		cellNum = (cellNum + 1) % this.cells.length
		if(cellNum == startingCellNumber) {
			// we have gone through the whole list with no valid neighbors
			return false
		}
		
	}
}

CellGroup.prototype.getTopLeft = function () {
	var cells = this.cells, topLeftCell = cells[0]
	
	for (var i = 1; i < cells.length; i++) {
		if (cells[i].x <= topLeftCell.x) {
			if (cells[i].y < topLeftCell.y) topLeftCell = cells[i]
		}
	}
	
	return topLeftCell
}

// A class for a single cell in the Kenken which houses data on the cell and methods for finding adjacent cells
function Cell (kenken, value) {
    this.kenken = kenken
    this.x = undefined
    this.y = undefined
	this.cellGroup = undefined
	this.value = value
}

// Function for setting the cell group that a cell belongs to
Cell.prototype.setCellGroup = function(cellGroup) {
	this.cellGroup = cellGroup
}

// Return an array of the cells neighbors
Cell.prototype.getNeighbors = function () {
    var neighbors = []
    
    if (this.kenken.settings.torus) {
    	if (this.x > 0) neighbors.push(this.kenken.board[this.x-1][this.y])
    	else neighbors.push(this.kenken.board[this.kenken.size-1][this.y])
	    if (this.y > 0) neighbors.push(this.kenken.board[this.x][this.y-1])
    	else neighbors.push(this.kenken.board[this.x][this.kenken.size-1])
	    if (this.x < this.kenken.size - 1) neighbors.push(this.kenken.board[this.x+1][this.y])
    	else neighbors.push(this.kenken.board[0][this.y])
	    if (this.y < this.kenken.size - 1) neighbors.push(this.kenken.board[this.x][this.y+1])
    	else neighbors.push(this.kenken.board[this.x][0])
    }
    else
    {
	    if (this.x > 0) neighbors.push(this.kenken.board[this.x-1][this.y])
	    if (this.y > 0) neighbors.push(this.kenken.board[this.x][this.y-1])
	    if (this.x < this.kenken.size - 1) neighbors.push(this.kenken.board[this.x+1][this.y])
	    if (this.y < this.kenken.size - 1) neighbors.push(this.kenken.board[this.x][this.y+1])
    }
    
    return neighbors
}

// Return an object with the cell's neighbors indexed by relative location
Cell.prototype.getNeighborsOriented = function () {
    var neighbors = {}
    
    if (this.kenken.settings.torus) {
    	neighbors.left=this.x > 0 ? this.kenken.board[this.x-1][this.y] : this.kenken.board[this.kenken.size-1][this.y]
    	neighbors.up = this.y > 0 ? this.kenken.board[this.x][this.y-1] : this.kenken.board[this.x][this.kenken.size-1]
	    neighbors.right = this.x < this.kenken.size - 1 ? this.kenken.board[this.x+1][this.y] : this.kenken.board[0][this.y]
		neighbors.down = this.y < this.kenken.size - 1 ? this.kenken.board[this.x][this.y+1] : this.kenken.board[this.x][0]
    }
    else
    {
	    if (this.x > 0) neighbors.left = this.kenken.board[this.x-1][this.y]
	    if (this.y > 0) neighbors.up = this.kenken.board[this.x][this.y-1]
	    if (this.x < this.kenken.size - 1) neighbors.right = this.kenken.board[this.x+1][this.y]
	    if (this.y < this.kenken.size - 1) neighbors.down = this.kenken.board[this.x][this.y+1]
    }
    
    return neighbors
}

// Function to shuffle an array, using to shuffle operations array.
function shuffleInputArray(array, seed) {
	// Randomly shuffle the array, doing the algorithm once forward and once
	// backward, to help create more "randomness"
	for (var i = 0; i < array.length-1; i++) {
		// Generate a random integer in the range [i,n-1]
		// Since seed.random() generates a number in the range [0,1)
		var randomNum = Math.floor((array.length-i)*seed.random()+i)
		
		//swap the array at spots i and randomNum
		var numToSwap = array[i]
		array[i] = array[randomNum]
		array[randomNum] = numToSwap
    }
}

// Function to generate an array with the numbers 1 through n in a random order
function shuffledNumberArray (n, seed) {
	var numberArray=[]
	// Fill the array with numbers 1 through n
	for(var i = 0; i < n; i++) {
		numberArray.push(i+1)
	}
	
	shuffleInputArray(numberArray,seed)
	
	//return then shuffled array
	return numberArray
}

/*
 * The following code is classes for the operations. Not sure if inheritance/abstract functions are a things in 
 * javascript so currently just making a different class for each of them.
 */
 
// The class for no operation, used for a single cell
function SingleCell() {
	this.minCells = 1 // Can operate on a minimum of 1 cell
	this.maxCells = 1 // Can operate on a maximum of 1 cell
	this.symbol = ''
}
 
// The operation function for a single cell
SingleCell.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}

	return arrayOfNumbers[0]
}
 
// The class for addition
function Addition() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = undefined // No max number of cells it can operate on
	this.symbol = '+'
}
 
// The operation function for addition
Addition.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	 
	var resultOfOperation = 0;
	for(var i = 0; i < arrayOfNumbers.length; i++) {
		resultOfOperation = arrayOfNumbers[i] + resultOfOperation
	}
	return resultOfOperation
}
 
// The class for subtraction
function Subtraction() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 2 // Can operate on a maximum of 2 cells
	this.symbol = '-'
}
 
// The operation function for subtraction
Subtraction.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	 
	// Should only have two numbers, subtract the smaller one from the larger one
	var resultOfOperation = false;
	if(arrayOfNumbers[0] > arrayOfNumbers[1]) {
		resultOfOperation = arrayOfNumbers[0] - arrayOfNumbers[1]
	} else {
		resultOfOperation = arrayOfNumbers[1] - arrayOfNumbers[0]
	}
	return resultOfOperation
 }
 
// The class for multiplication
function Multiplication() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = undefined // No max number of cells it can operate on
	this.symbol = '&times;'
 }
 
// The operation function for multiplication
Multiplication.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var resultOfOperation = 1;
	for(var i = 0; i < arrayOfNumbers.length; i++) {
		resultOfOperation = arrayOfNumbers[i] * resultOfOperation
	}
	return resultOfOperation
}
 
// The class for division
function Division() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 2 // Can operate on a maximum of 2 cells
	this.symbol = '&divide;'
}
 
// The operation function for division
Division.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	// Should only have two numbers, check to see if they divide evenly. If not, return 0 indicating division failure
	var resultOfOperation = false;
	if((arrayOfNumbers[0]/arrayOfNumbers[1])%1==0) {
		// If true, then this was an integer
		resultOfOperation = arrayOfNumbers[0]/arrayOfNumbers[1]
	} else if((arrayOfNumbers[1]/arrayOfNumbers[0])%1==0) {
		// If true, then this was an integer
		resultOfOperation = arrayOfNumbers[1]/arrayOfNumbers[0]
	}
	return resultOfOperation
}

// The class for min
function Minimum() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 3 // Can operate on a maximum of 3 cells
	this.symbol = 'min'
}

// The operation function for minimum
Minimum.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var resultOfOperation = arrayOfNumbers[0]
	for(var i = 1; i < arrayOfNumbers.length; i++) {
		if(resultOfOperation > arrayOfNumbers[i]) {
			resultOfOperation = arrayOfNumbers[i]
		}
	}
	return resultOfOperation
}

// The class for max
function Maximum() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 3 // Can operate on a maximum of 3 cells
	this.symbol = 'max'
}

// The operation function for maximum
Maximum.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var resultOfOperation = arrayOfNumbers[0]
	for(var i = 1; i < arrayOfNumbers.length; i++) {
		if(resultOfOperation < arrayOfNumbers[i]) {
			resultOfOperation = arrayOfNumbers[i]
		}
	}
	return resultOfOperation
}

// The class for range
function Range() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 4 // Can operate on a maximum of 4 cells
	this.symbol = 'range'
}

// The operation function for range
Range.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var min = arrayOfNumbers[0]
	var max = arrayOfNumbers[0]
	for(var i = 1; i < arrayOfNumbers.length; i++) {
		if(max < arrayOfNumbers[i]) {
			max = arrayOfNumbers[i]
		} else if(min > arrayOfNumbers[i]) {
			min = arrayOfNumbers[i]
		}
	}
	return max-min
}

// The class for mod
function Mod() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 2 // Can operate on a maximum of 2 cells
	this.symbol = '%'
}

// The operation function for mod
Mod.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	// Should only have two numbers, take larger % smaller
	var resultOfOperation = false;
	if(arrayOfNumbers[0] > arrayOfNumbers[1]) {
		resultOfOperation = arrayOfNumbers[0] % arrayOfNumbers[1]
	} else {
		resultOfOperation = arrayOfNumbers[1] % arrayOfNumbers[0]
	}
	return resultOfOperation
}

// The class for average
function Average() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 4 // Can operate on a maximum of 4 cells
	this.symbol = 'avg'
}

// The operation function for average
Average.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var sum = 0;
	//Calculate the average of the values
	for(var i = 0; i < arrayOfNumbers.length; i++) {
		sum = sum + arrayOfNumbers[i]
	}
	var average = sum / arrayOfNumbers.length
	if((average % 1)==0) {
		return average
	} else {
		return false
	}
}

// The class for parity
function Parity() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 2 // Can operate on a maximum of 2 cells
	this.symbol = 'par'
}

// The operation function for parity
Parity.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	var firstNum = arrayOfNumbers[0]
	var secondNum = arrayOfNumbers[1]
	if( (firstNum + secondNum) % 2 == 0) {
		// They have the same parity, return 1
		return 1
	} else {
		// They have different parity, return 0
		return 0
	}
}

// The class for gcd
function Gcd() {
	this.minCells = 2 // Can operate on a minimum of 2 cells
	this.maxCells = 3 // Can operate on a maximum of 4 cells
	this.symbol = 'gcd'
}

// The operation function for gcd
Gcd.prototype.operation = function(arrayOfNumbers) {
	// Verify that array length is within the size constraints, if not return 0 meaning the operation failed
	if(arrayOfNumbers.length > this.maxCells || arrayOfNumbers.length < this.minCells) {
		return false
	}
	
	// A recursive function for finding the gcd of two numbers
	function gcd(a,b) { return b ? gcd(b,a%b) : a }
	
	var currentGCD = arrayOfNumbers[0]
	for(var i = 1; i < arrayOfNumbers.length; i++) {
		currentGCD = gcd(currentGCD,arrayOfNumbers[i])
	}
	return currentGCD
}

var optionKeys = ["size", "difficulty", "seed", "maxGroupSize", "torus"], optionOperations = ["addition", "subtraction", "multiplication", "division", "min", "max", "range", "mod", "gcd", "par", "avg"]

function encodeOptions (options) {
    var data = [], operations = []
    
    for (var i = 0; i < optionKeys.length; i++) {
        data.push(options[optionKeys[i]])
    }
    
    for (var i = 0; i < optionOperations.length; i++) {
        operations.push(options.operations[optionOperations[i]]? 1 : 0)
    }
    
    console.log(options, data, operations)
    return btoa(JSON.stringify(data) + "::" + JSON.stringify(operations))
}

function decodeOptions (string) {
    string = atob(string).split("::")
    var data = JSON.parse(string[0]), operations = JSON.parse(string[1]), options = {operations: {}}
    
    for (var i = 0; i < optionKeys.length; i++) {
        options[optionKeys[i]] = data[i]
    }
    
    for (var i = 0; i < optionOperations.length; i++) {
        options.operations[optionOperations[i]] = operations[i] ? true : false
    }
    
    return options
}


currentKenken = generateKenken({
    size: 5,
    operations: {
        addition: true,
        subtraction: false,
        multiplication: false,
        division: false,
        min: false,
        max: false,
        range: false,
        mod: false,
        avg: false,
        par: false,
        gcd: false
    },
    difficulty: 2,
    maxGroupSize: 3,
    torus: false
})


// currentKenken = generateKenken({
//     size: $("#size").spinner("value"),
//     operations: {
//         addition: $("#addition").is(":checked"),
//         subtraction: $("#subtraction").is(":checked"),
//         multiplication: $("#multiplication").is(":checked"),
//         division: $("#division").is(":checked"),
//         min: $("#min").is(":checked"),
//         max: $("#max").is(":checked"),
//         range: $("#range").is(":checked"),
//         mod: $("#mod").is(":checked"),
//         avg: $("#avg").is(":checked"),
//         par: $("#par").is(":checked"),
//         gcd: $("#gcd").is(":checked")
//     },
//     difficulty: $("#difficulty").spinner("value"),
//     maxGroupSize: $("#groupSize").spinner("value"),
//     torus: $("#torus").is(":checked")
// })


console.log(currentKenken)

function print_cell_groups(kenken){
    var size = kenken.size

    for (var y = 0; y < size; y++){
        for (var x = 0; x < size; x++) {
            process.stdout.write(kenken.board[x][y].cellGroup.groupID + " ")
        }
        console.log("\n")
    }
} 

function get_cell_groups(kenken){
    var size = kenken.size
    out = []
    for (var y = 0; y < size; y++){
        row = []
        for (var x = 0; x < size; x++) {
            row.push(kenken.board[x][y].cellGroup.groupID)
        }
        out.push(row)
    }
    return out
} 

function print_sol(kenken){
    var size = kenken.size

    for (var y = 0; y < size; y++){
        for (var x = 0; x < size; x++) {
            process.stdout.write(kenken.board[x][y].value + " ")
        }
        console.log("\n")
    }
} 

function get_sol(kenken){
    var size = kenken.size
    const out = []
    for (var y = 0; y < size; y++){
        row = []
        for (var x = 0; x < size; x++) {
            row.push(kenken.board[x][y].value)
        }
        out.push(row)
    }
    return out
} 

function print_cell_group_info(kenken){
    var size = kenken.cellGroups.length

    for (var y = 0; y < size; y++){
     
        console.log(kenken.cellGroups[y].groupID + "--" + kenken.cellGroups[y].operationDescription)
    }
}


function get_cell_group_info(kenken){

    const cell_groups = []

    var size = kenken.cellGroups.length

    for (var y = 0; y < size; y++){
        cell_groups.push([kenken.cellGroups[y].groupID,kenken.cellGroups[y].operationDescription])
    }

    return cell_groups
}

// print_sol(currentKenken)
// print_cell_groups(currentKenken)
// print_cell_group_info(currentKenken)

// console.log(get_cell_group_info(currentKenken))
// console.log(get_sol(currentKenken))
// console.log(get_cell_groups(currentKenken))

// ken_json = {
//     'solution' : get_sol(currentKenken),
//     'cellGroups': get_cell_groups(currentKenken),
//     'cellValues': get_cell_group_info(currentKenken)
// }

// console.log(JSON.stringify(ken_json))

const kenken_set = new Set()
const num_games = 10000000
const games_json = []
for (var i = 0; i< num_games; i++){
    currentKenken = generateKenken({
        size: 3,
        operations: {
            addition: true,
            subtraction: false,
            multiplication: false,
            division: false,
            min: false,
            max: false,
            range: false,
            mod: false,
            avg: false,
            par: false,
            gcd: false
        },
        difficulty: 1,
        maxGroupSize: 3,
        torus: false
    })

    ken_json = {
        'solution' : get_sol(currentKenken),
        'cellGroups': get_cell_groups(currentKenken),
        'groupValues': get_cell_group_info(currentKenken)
    }
    
    str = JSON.stringify(ken_json)

    if (!(kenken_set.has(str))){
        kenken_set.add(str)
        games_json.push(ken_json)
    }
}

console.log("The num unique games generated: ")
console.log(games_json.length)
let data = JSON.stringify(games_json);
let file_name = "111games_json-4x4_cage2.json"
console.log("Saving the games json file here: " + file_name)
fs.writeFileSync(file_name, data);