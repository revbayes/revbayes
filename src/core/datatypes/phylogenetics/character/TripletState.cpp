//#include "TripletState.h"
//#include "DnaState.h"
//#include "RbException.h"
//#include "StringUtilities.h"
//
//#include <sstream>
//
//using namespace RevBayesCore;
//
///** Default constructor */
//TripletState::TripletState(size_t n) : DiscreteCharacterState( 64 ),
//    is_gap( false ),
//    is_missing( false ),
//    index_single_state( 0 ),
//    num_observed_states( 0 ),
//    state(64)
//{
//    
//    setState("---");
//}
//
//
///** Constructor that sets the observation */
//TripletState::TripletState(const std::string &s) : DiscreteCharacterState( 64 ),
//    is_gap( false ),
//    is_missing( false ),
//    index_single_state( 0 ),
//    num_observed_states( 0 ),
//    state(64)
//{
//    
//    setState(s);
//}
//
//
///* Clone object */
//TripletState* TripletState::clone(void) const
//{
//    
//    return new TripletState( *this );
//}
//
//AminoAcidState TripletState::getAminoAcidState( void ) const
//{
//    
//    AminoAcidState aa;
//    
//    std::vector<unsigned int> Triplet_pos = getTripletStates();
//    
//    // test if the first Triplet position is an 'A'
//    if ( Triplet_pos[0] == 0 )
//    {
//        
//        // test if the second Triplet position is an 'A'
//        if ( Triplet_pos[1] == 0 )
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a Lysine
//                aa.addState("K");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Asparagine
//                aa.addState("N");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 1 ) // test if the second Triplet position is an 'C'
//        {
//            
//            // we have a Threonine
//            aa.addState("T");
//            
//        }
//        else if ( Triplet_pos[1] == 2 ) // test if the second Triplet position is an 'G'
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a Arginine
//                aa.addState("R");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Serine
//                aa.addState("S");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 3 ) // test if the second Triplet position is an 'T'
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'C' or a 'T'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Isoleucine
//                aa.addState("I");
//            }
//            else if ( Triplet_pos[2] == 2 )
//            {
//                // we have a Methionine
//                aa.addState("M");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else
//        {
//            throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//        }
//        
//    }
//    else if ( Triplet_pos[0] == 1 ) // test if the first Triplet position is a 'C'
//    {
//        
//        // test if the second Triplet position is an 'A'
//        if ( Triplet_pos[1] == 0 )
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a Glutamine
//                aa.addState("Q");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Histidine
//                aa.addState("H");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 1 ) // test if the second Triplet position is an 'C'
//        {
//            
//            // we have a Proline
//            aa.addState("P");
//            
//        }
//        else if ( Triplet_pos[1] == 2 ) // test if the second Triplet position is an 'G'
//        {
//            
//            // we have a Arginine
//            aa.addState("R");
//            
//        }
//        else if ( Triplet_pos[1] == 3 ) // test if the second Triplet position is an 'T'
//        {
//            
//            // we have a Leucine
//            aa.addState("L");
//            
//        }
//        else
//        {
//            throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//        }
//        
//    }
//    else if ( Triplet_pos[0] == 2 ) // test if the first Triplet position is a 'G'
//    {
//        
//        // test if the second Triplet position is an 'A'
//        if ( Triplet_pos[1] == 0 )
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2]== 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a Glutamic Acid
//                aa.addState("E");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Aspartic Acid
//                aa.addState("D");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 1 ) // test if the second Triplet position is an 'C'
//        {
//            
//            // we have a Alanine
//            aa.addState("A");
//            
//        }
//        else if ( Triplet_pos[1] == 2 ) // test if the second Triplet position is an 'G'
//        {
//            
//            // we have a Glycine
//            aa.addState("G");
//            
//        }
//        else if ( Triplet_pos[1] == 3 ) // test if the second Triplet position is an 'T'
//        {
//            
//            // we have a Valine
//            aa.addState("V");
//            
//        }
//        else
//        {
//            throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//        }
//        
//    }
//    else if ( Triplet_pos[0] == 3 ) // test if the first Triplet position is an 'T'
//    {
//        
//        // test if the second Triplet position is an 'A'
//        if ( Triplet_pos[1] == 0 )
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a stopp Triplet
//                aa.setMissingState( true );
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Tyrosine
//                aa.addState("Y");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 1 ) // test if the second Triplet position is an 'C'
//        {
//            
//            // we have a Serine
//            aa.addState("S");
//            
//        }
//        else if ( Triplet_pos[1] == 2 ) // test if the second Triplet position is an 'G'
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 )
//            {
//                // we have a stopp Triplet
//                aa.setMissingState( true );
//            }
//            else if ( Triplet_pos[2] == 2 )
//            {
//                // we have a Tryptophan
//                aa.addState("W");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Cysteine
//                aa.addState("C");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 3 ) // test if the second Triplet position is an 'T'
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a Leucine
//                aa.addState("L");
//            }
//            else if ( Triplet_pos[2] == 1 || Triplet_pos[2] == 3 )
//            {
//                // we have a Phenylalanine
//                aa.addState("F");
//            }
//            else
//            {
//                throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//            }
//            
//        }
//        else
//        {
//            throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//        }
//        
//    }
//    else
//    {
//        throw RbException() << "Unknown amino acid for Triplet '" <<  getStringValue() << "'";
//    }
//    
//    
//    return aa;
//}
//
//
//std::string TripletState::getDataType( void ) const
//{
//    return "Triplet";
//}
//
//
//std::string TripletState::getStateLabels( void ) const
//{
//    
//    static const std::string stateLabels = "AAA...TTT";
//    
//    return stateLabels;
//}
//
//std::string TripletState::getStringValue(void) const
//{
//    
//    if ( isMissingState() )
//    {
//        return "?";
//    }
//    
//    if ( isGapState() )
//    {
//        return "-";
//    }
//    
//    std::vector<unsigned int> Triplet_pos = getTripletStates();
//    
//    std::string str_val = "";
//    
//    for (size_t i=0; i<3; ++i)
//    {
//        
//        switch ( Triplet_pos[i] )
//        {
//            case 0x0:
//                str_val += "A";
//                break;
//            case 0x1:
//                str_val += "C";
//                break;
//                //            case 0x3:
//                //                str_val += "M";
//                //                break;
//            case 0x2:
//                str_val += "G";
//                break;
//                //            case 0x5:
//                //                str_val += "R";
//                //                break;
//                //            case 0x6:
//                //                str_val += "S";
//                //                break;
//                //            case 0x7:
//                //                str_val += "V";
//                //                break;
//            case 0x3:
//                str_val += "T";
//                break;
//                //            case 0x9:
//                //                str_val += "W";
//                //                break;
//                //            case 0xA:
//                //                str_val += "Y";
//                //                break;
//                //            case 0xB:
//                //                str_val += "H";
//                //                break;
//                //            case 0xC:
//                //                str_val += "K";
//                //                break;
//                //            case 0xD:
//                //                str_val += "D";
//                //                break;
//                //            case 0xE:
//                //                str_val += "B";
//                //                break;
//            case 0xF:
//                str_val += "N";
//                break;
//                //            case 0xF:
//                //                str_val += "-";
//                //                break;
//                
//            default:
//                str_val += "?";
//        }
//    }
//    
//    return str_val;
//}
//
//
//std::vector<unsigned int> TripletState::getTripletStates( void ) const
//{
//    
//    
//    std::vector<unsigned int> Triplet_pos = std::vector<unsigned int>(3,0);
//    Triplet_pos[0] = getStateIndex() % 4;
//    Triplet_pos[1] = int(getStateIndex() / 4) % 4;
//    Triplet_pos[2] = int(getStateIndex() / 16) % 4;
//    
//    return Triplet_pos;
//}
//
//
//bool TripletState::isStopTriplet( void ) const
//{
//    
//    std::vector<unsigned int> Triplet_pos = getTripletStates();
//    
//    
//    if ( Triplet_pos[0] == 3 ) // test if the first Triplet position is an 'T'
//    {
//        
//        // test if the second Triplet position is an 'A'
//        if ( Triplet_pos[1] == 0 )
//        {
//            
//            // test if the third Triplet position is an 'A' or a 'G'
//            if ( Triplet_pos[2] == 0 || Triplet_pos[2] == 2 )
//            {
//                // we have a stopp Triplet
//                return true;
//            }
//            
//        }
//        else if ( Triplet_pos[1] == 2 ) // test if the second Triplet position is an 'G'
//        {
//            
//            // test if the third Triplet position is an 'A'
//            if ( Triplet_pos[2] == 0 )
//            {
//                // we have a stopp Triplet
//                return true;
//            }
//            
//        }
//        
//    }
//    
//    
//    return false;
//}
//
//
//bool TripletState::isGapState( void ) const
//{
//    return is_gap;
//}
//
//
//bool TripletState::isMissingState( void ) const
//{
//    return is_missing;
//}
//
//
//void TripletState::setGapState( bool tf )
//{
//    is_gap = tf;
//}
//
//
//void TripletState::setMissingState( bool tf )
//{
//    is_missing = tf;
//}
//
//
//void TripletState::setState(const std::string &s)
//{
//    std::string symbol = s;
//    if ( symbol == "???" )
//    {
//        setMissingState( true );
//    }
//    else if ( symbol == "---" )
//    {
//        setGapState( true );
//    }
//    else
//    {
//        
//        /* A C G T */
//        StringUtilities::replaceSubstring(symbol, "U", "T");
//        DnaState dna_pos_0 = DnaState( std::string(1, symbol[0]) );
//        DnaState dna_pos_1 = DnaState( std::string(1, symbol[1]) );
//        DnaState dna_pos_2 = DnaState( std::string(1, symbol[2]) );
//        
//        RbBitSet bs_pos_0 = dna_pos_0.getState();
//        RbBitSet bs_pos_1 = dna_pos_1.getState();
//        RbBitSet bs_pos_2 = dna_pos_2.getState();
//        
//        
//        num_observed_states = 0;
//        state.reset();
//        
//        for (size_t i=0; i<4; ++i)
//        {
//            
//            // test if the bit is set for the first Triplet position
//            if ( bs_pos_0.test( i ) == true )
//            {
//                
//                for (size_t j=0; j<4; ++j)
//                {
//                    
//                    // test if the bit is set for the second Triplet position
//                    if ( bs_pos_1.test( j ) == true )
//                    {
//                        
//                        for (size_t k=0; k<4; ++k)
//                        {
//                            
//                            // test if the bit is set for the third Triplet position
//                            if ( bs_pos_2.test( k ) == true )
//                            {
//                                
//                                ++num_observed_states;
//                                size_t current_state = i*16 + j*4 + k;
//                                state.set( current_state );
//                                index_single_state = current_state;
//                                
//                            }
//                            
//                        } // end for-loop over all possible states for the third Triplet position
//                        
//                    }
//                    
//                } // end for-loop over all possible states for the second Triplet position
//                
//            }
//            
//        } // end for-loop over all possible states for the first Triplet position
//        
//        
//    } // end if this is not a missing or gap state
//    
//}
//
//
//void TripletState::addState(const std::string &symbol)
//{
//    ++num_observed_states;
//    
//    std::string labels = getStateLabels();
//    size_t pos = labels.find(symbol);
//    
//    state.set( pos );
//    index_single_state = pos;
//}
//
//
//
//RbBitSet TripletState::getState(void) const
//{
//    return state;
//}
//
//
//void TripletState::setToFirstState(void)
//{
//    num_observed_states = 1;
//    index_single_state = 0;
//    state.reset();
//    state.set( 0 );
//}
//
//
//void TripletState::setStateByIndex(size_t index)
//{
//    
//    num_observed_states = 1;
//    index_single_state = index;
//    state.reset();
//    state.set( index );
//}
//
//
