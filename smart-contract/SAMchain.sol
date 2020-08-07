pragma solidity ^0.5.8;
pragma experimental ABIEncoderV2;

contract SAMchain {

// 15 June 2020
// GAMZE GURSOY & CHARLOTTE BRANNON
// YALE UNIVERSITY
//**********************************************************************************************
// DEFINE STRUCTS
//**********************************************************************************************

    // This structure is how the data should be returned from the query function.
    // You do not have to store relations this way in your contract, only return them.
    // geneName and drugName must be in the same capitalization as it was entered. E.g. if the original entry was GyNx3 then GYNX3 would be considered incorrect.
    // Percentage values must be acurrate to 6 decimal places and will not include a % sign. E.g. "35.123456"
    struct ReadStruct {
        string qname;
        uint flag;
        string rname;
        uint position;
        uint mapq;
        string cigar;
        uint rnext;
        uint pnext;
        uint tlen;
        string seq;
        string qual;
        uint index;
    }

//**********************************************************************************************
// INITIALIZE MAPPINGS AND ARRAYS
//**********************************************************************************************

    // mappings and corresponding arrays to hold data values and their keys, respectively
    mapping (string => uint[]) qnames; //key = geneName
    mapping (uint => ReadStruct) database;
    mapping (uint => uint[]) positions;
    mapping (string => mapping(uint => uint[])) chrPositions;

    // counter to assign each entry a unique index
    uint rcounter = 0;
    uint ecounter = 1;
    // parameters for toString()
    uint prec_all = 8; 
    uint prec_digits = 2;


//**********************************************************************************************
// CORE FUNCTIONS
//**********************************************************************************************

    /*  Insert an observation into your contract, following the format defined in the data readme. 
        This function has no return value. If it completes it will be assumed the observations was recorded successfully. 
        Note: case matters for geneName and drugName. GyNx3 and gynx3 are treated as different genes.
    */
    function insertRead (
        string memory qname,
        uint flag,
        string memory rname,
        uint position,
        uint mapq,
        string memory cigar,
        uint rnext,
        uint pnext,
        uint tlen,
        string memory seq,
        string memory qual
    ) public { 

        uint index;
        index = rcounter;
        database[index] = ReadStruct(qname, flag, rname, position, mapq, cigar, rnext, pnext, tlen, seq, qual, index);
        qnames[qname].push(rcounter);
        positions[position].push(rcounter);
        chrPositions[rname][position].push(rcounter);
        rcounter++;
    }

//Retrieves entire read database
function returnDatabase() public view returns(ReadStruct[] memory){
        ReadStruct[] memory contents = new ReadStruct[](rcounter);
        for(uint i; i < rcounter; i++) {
           contents[i] = database[i];
        }
       return contents;
    }

// Retrieves reads mapping to a point POS value
function queryPosition(
        uint position
    ) public view returns (ReadStruct[] memory) {
        ReadStruct[] memory contents = new ReadStruct[](rcounter);
        uint[] memory positionSearch;
        ReadStruct[] memory empty;
        // if database is empty, return empty array
        if (rcounter == 0) {
            return empty;
        }
        positionSearch = positions[position];
        // initialize memory structs and variables
        for (uint i; i < positionSearch.length; i++) {
            contents[i] = database[positionSearch[i]];
        }
        return contents;
    }

// Retrieves reads mapping to a point POS value
function queryChrPosition(
        string memory chr,
        uint position
    ) public view returns (ReadStruct[] memory) {
        ReadStruct[] memory contents = new ReadStruct[](rcounter);
        uint[] memory chrPositionSearch;
        ReadStruct[] memory empty;
        // if database is empty, return empty array
        if (rcounter == 0) {
            return empty;
        }
        chrPositionSearch = chrPositions[chr][position];
        // initialize memory structs and variables
        for (uint i; i < chrPositionSearch.length; i++) {
            contents[i] = database[chrPositionSearch[i]];
        }
        return contents;
    }

// Retrieves reads mapping to a point POS value
function queryQname(
        string memory qname
    ) public view returns (ReadStruct[] memory) {
        ReadStruct[] memory contents = new ReadStruct[](rcounter);
        uint[] memory qnameSearch;
        ReadStruct[] memory empty;
        // if database is empty, return empty array
        if (rcounter == 0) {
            return empty;
        }
        qnameSearch = qnames[qname];
        // initialize memory structs and variables
        for (uint i; i < qnameSearch.length; i++) {
            contents[i] = database[qnameSearch[i]];
        }
        return contents;
    }
//**********************************************************************************************
// BASIC UTILITIES
//**********************************************************************************************

    /*  CMB: function to convert uints to strings
        from https://github.com/willitscale/solidity-util/blob/master/lib/Integers.sol
        edited by CMB for the percentage string case.
    */
    function toString(
        uint _base, 
        uint numDigits, 
        bool lessThanOne, 
        uint precision
        ) internal pure returns (string memory) {
        bytes memory _tmp = new bytes(32);
        uint i;
        for(i; _base > 0; i++) {
            _tmp[i] = byte(uint8((_base % 10) + 48));
            _base /= 10;
        }
        bytes memory _real = new bytes(i--);
        for(uint j; j < _real.length; j++) {
            _real[j] = _tmp[i--];
        }
        bytes memory deci_real = new bytes(_real.length + numDigits);
        if (precision != 0) {
            uint count;
            uint tally;
            for (uint k; k < deci_real.length; k++) {
                if (lessThanOne == true) {
                    if (k < precision) {
                        if (k == 0) {
                            deci_real[count] = "0";
                            count++;
                        } else if (k ==1) {
                            deci_real[count] = ".";
                            count++;
                        } else if (k!=0 && k!= 1 && k < (precision - (numDigits - 1))) {
                            deci_real[count] = "0";
                            count++;
                        } else if (k!=0 && k!= 1 && k >= (precision - (numDigits - 1))) {
                            deci_real[count] = _real[tally];
                            tally++;
                            count++;
                        }
                    }
                } else {
                    if (k == numDigits) {
                        deci_real[count] = ".";
                        count++;
                        deci_real[count] = _real[k];
                        count++;
                    } else if (k < _real.length) {
                        deci_real[count] = _real[k];
                        count++;
                    }
                }
            }
            bytes memory bytesContainerTrimmed = new bytes(count);
            for (uint256 charCounter; charCounter < count; charCounter++) {
                bytesContainerTrimmed[charCounter] = deci_real[charCounter];
            }
            return string(bytesContainerTrimmed);
        }
        return string(_real);
    }

    /*  CMB: function to convert uints to strings
        from https://github.com/willitscale/solidity-util/blob/master/lib/Integers.sol
    */
    function toStringSimple(uint _base) internal pure returns (string memory) {
        bytes memory _tmp = new bytes(32);
        uint i;
        for(i; _base > 0; i++) {
            _tmp[i] = byte(uint8((_base % 10) + 48));
            _base /= 10;
        }
        bytes memory _real = new bytes(i--);
        for(uint j; j < _real.length; j++) {
            _real[j] = _tmp[i--];
        }
        return string(_real);
    }

    /*  CMB: function to return the number of digits in a uint
        from https://github.com/willitscale/solidity-util/blob/master/lib/Integers.sol
    */
    function numDigits(uint number) internal pure returns (uint result) {
        uint digits;
        while (number != 0) {
            number /= 10;
            digits++;
        }
        return digits;
    }

    /*  CMB: function to divide two uints and return an answer with high precision
    */
    function percent(
        uint numerator, 
        uint denominator, 
        uint precision
        ) internal pure returns(uint quotient) {
        uint _numerator  = numerator * 10 ** (precision+1);
        uint _quotient =  ((_numerator / denominator + 5) / 10);
        return (_quotient);
    }

    /*  CMB: function to convert a string to a uint
        from https://github.com/willitscale/solidity-util/blob/master/lib/Integers.sol
    */
    function toInt(string memory _value) internal pure returns (uint _ret) {
        bytes memory _bytesValue = bytes(_value);
        uint j = 1;
        for(uint i = _bytesValue.length-1; i >= 0 && i < _bytesValue.length; i--) {
            assert(uint8(_bytesValue[i]) >= 48 && uint8(_bytesValue[i]) <= 57);
            _ret += (uint8(_bytesValue[i]) - 48)*j;
            j*=10;
        }
    }

    /*  CMB: function to compare two strings
    */
    function compareStrings(
        string memory a, 
        string memory b
        ) internal pure returns (bool){ 
        return (keccak256(abi.encodePacked((a))) == keccak256(abi.encodePacked((b)))); 
    }

    /*  CMB: function to convert a string to bytes32 data type
    */
    function toBytes32(string memory _string) internal pure returns (bytes32) {
        bytes32 _stringBytes;
        assembly {
        _stringBytes := mload(add(_string, 32))
        }
        return _stringBytes;
    }
    
    /*  CMB: function to convert a bytes32 to string data type
    */
    function toShortString32(bytes32 _data) internal pure returns (string memory) {
        bytes memory _bytesContainer = new bytes(32);
        uint256 _charCount;
        for (uint256 _bytesCounter; _bytesCounter < 32; _bytesCounter++) {
            bytes1 _char = bytes1(bytes32(uint256(_data) * 2 ** (8 * _bytesCounter)));
            if (_char != 0) {
                _bytesContainer[_charCount] = _char;
                _charCount++;
            }
        }
        bytes memory _bytesContainerTrimmed = new bytes(_charCount);
        for (uint256 _charCounter; _charCounter < _charCount; _charCounter++) {
            _bytesContainerTrimmed[_charCounter] = _bytesContainer[_charCounter];
        }
        return string(_bytesContainerTrimmed);
    }

} // END OF CONTRACT