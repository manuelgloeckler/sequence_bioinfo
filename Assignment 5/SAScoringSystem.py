
__author__ = "Yanpeng Li, Josua Stadelmaier, Finn Mier"

import logging

logger = logging.getLogger()

class ScoringSystem:
    '''Abstract class describing a scoring system. None of its methods are actually implemented.
    '''
    def compare(self, a, b) -> int:
        ''' Returns a score dependening on the similarity of a and b.

        Parameters:
            a : first value to be compared. Type must be supported by the scoring method.
            b : second value to be compared. Type must be supported by the scoring method.
        Returns:
            Score resulting from the comparison if the two given values.
        '''
        raise NotImplementedError
    
    def gap_penalty(self, length : int) -> int:
        ''' Returns a score depending on the lenght of the gap.

        Parameters:
            length (int): length of the gap.
        Returns:
            Score computed for the gap with the given length.
        '''
        raise NotImplementedError

class EqualLinearScoring(ScoringSystem):
    ''' Class implementing a scoring system which punishes any mismatch equally, any
    match equally and has a linear gap cost.
        
    Paramters:
        match_score (int): reward for aligning two equal values.
        mismatch_score (int): "reward" for aliging two different values.
        gap_score (int): score which is linearly scaled with gap length.
    '''
    def __init__(self, match_score : int, mismatch_score : int, gap_score : int):
        # sanity check
        self._sanity_check(match_score, mismatch_score, gap_score)
        self.match_score = match_score
        self.mismatch_score = mismatch_score    
        self.gap_score = gap_score
        logger.info("Using equal scoring with linear gap function and following scores: match: {}, mismatch: {}, gap penalty: {}".format(match_score, mismatch_score, gap_score))
    
    def _sanity_check(self, match_score : int, mismatch_score : int, gap_score : int):
        if match_score < mismatch_score:
            logger.error("Match score: {} was lower than mismatch score: {}".format(match_score, mismatch_score))
            raise ValueError
        if match_score < -gap_score:
            logger.error("Match score: {} was lower than gap score: {}".format(match_score, -gap_score))
            raise ValueError


    def compare(self, a, b) -> int:
        '''Compares the values a and b and returns this scoring systems match_score if they were equal
        and its mismatch_score if they weren't.

        Parameters:
            a : first value to be compared.
            b : second value to be compared.
        Returns:
            Score resulting from the comparison if the two given values.
        ''' 
        return self.match_score if a == b else self.mismatch_score

    def gap_penalty(self, length : int) -> int:
        '''Takes the length of a gap and returns the score penalty for increasing it by one. In this case the score is computed by
        this scoring systems gap_score times the given length.

        Parameters:
            length (int): length of the gap.
        Returns:
            Score computed for the gap with the given length.
        '''
        return -self.gap_score

class MatrixLinearScoring(ScoringSystem):
    ''' Class implementing a scoring system which looks up the scores for two different values in a scoring matrix and has a linear gap cost.
        
    Paramters:
        scoring_alphabet ({object : int}): Dictionary mapping values which are to be scored to integers.
        scoring_matrix_path ([[int]]): List of list containg scores to compare two values in scoring_alphabet.
        gap_score (int): score which is linearly scaled with gap length.
    '''
    def __init__(self, scoring_alphabet, scoring_matrix, gap_score : int):
        self._sanity_check(scoring_matrix, gap_score)
        self.scoring_alphabet = scoring_alphabet
        self.scoring_matrix = scoring_matrix    
        self.gap_score = gap_score
        logger.info("Using scoring matrix with linear gap function and following gap penalty: {}".format(gap_score))
    
    def _sanity_check(self, scoring_matrix, gap_score : int):
        lowest_score = min([min(l) for l in scoring_matrix])
        if lowest_score < -gap_score:
            logger.error("Lowest matching score: {} was lower than gap score: {}".format(lowest_score, -gap_score))
            raise ValueError

    def compare(self, a, b) -> int:
        '''Compares the values a and b and returns this scoring systems match_score if they were equal
        and its mismatch_score if they weren't.

        Parameters:
            a : first value to be compared.
            b : second value to be compared.
        Returns:
            Score resulting from the comparison if the two given values.
        ''' 
        return self.scoring_matrix[self.scoring_alphabet[a]][self.scoring_alphabet[b]]

    def gap_penalty(self, length : int) -> int:
        '''Takes the length of a gap and returns the score penalty for increasing it by one. In this case the score is computed by
        this scoring systems gap_score times the given length.

        Parameters:
            length (int): length of the gap.
        Returns:
            Score resulting from the comparison if the two given values.
        '''
        return -self.gap_score


class EqualAffineScoring(ScoringSystem):
    ''' Class implementing a scoring system which punishes any mismatch equally, any
    match equally and has an affine gap cost.
        
    Paramters:
        match_score (int): reward for aligning two equal values.
        mismatch_score (int): "reward" for aliging two different values.
        gap_open (int): penalty for opening a gap.
        gap_extend (int): penalty for  a gap.
    '''
    def __init__(self, match_score : int, mismatch_score : int, gap_open : int, gap_extend : int):
        self.match_score = match_score
        self.mismatch_score = mismatch_score    
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        logger.info("Using equal scoring with affine gap function and following scores: match: {}, mismatch: {}, gap open penalty: {}, gap extension penalty: {}".format(match_score, mismatch_score, gap_open, gap_extend))
     
    def _sanity_check(self, match_score : int, mismatch_score : int, gap_open : int, gap_extend : int):
        if match_score < mismatch_score:
            logger.error("Match score: {} was lower than mismatch score: {}".format(match_score, mismatch_score))
            raise ValueError
        
        if mismatch_score < - gap_open - gap_extend:
            logger.error("Mismatch score: {} was lower than gap open + gap extend: {}".format(mismatch_score, -gap_open - gap_extend))
            raise ValueError

        if gap_open < gap_extend:
            logger.error("Gap extend: {} was greater than gap open: {}.".format(gap_extend, gap_open))
            raise ValueError

    def compare(self, a, b) -> int:
        '''Compares the values a and b and returns this scoring systems match_score if they were equal
        and its mismatch_score if they weren't.

        Parameters:
            a : first value to be compared.
            b : second value to be compared.
        Returns:
            Score computed for increasing the gap with the given length.
        ''' 
        return self.match_score if a == b else self.mismatch_score

    def gap_penalty(self, length : int) -> int:
        '''Takes the length of a gap and returns the score penalty for increasing it by one. In this case the score is computed by
        this scoring systems gap_score times the given length.

        Parameters:
            length (int): length of the gap.
        Returns:
            Score computed for increasing the gap with the given length.
        '''
        if length == 1: return -self.gap_open
        elif length > 1: return -self.gap_extend

class MatrixAffineScoring(ScoringSystem):
    ''' Class implementing a scoring system which looks up the scores for two different values in a scoring matrix and has a affine gap cost.
        
    Paramters:
        scoring_alphabet ({object : int}): Dictionary mapping values which are to be scored to integers.
        scoring_matrix_path ([[int]]): List of list containg scores to compare two values in scoring_alphabet.
        gap_open (int): penalty for opening a gap.
        gap_extend (int): penalty for  a gap.
    '''
    def __init__(self, scoring_alphabet, scoring_matrix, gap_open : int, gap_extend : int):
        self._sanity_check(scoring_matrix, gap_open, gap_extend)
        self.scoring_alphabet = scoring_alphabet
        self.scoring_matrix = scoring_matrix    
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        logger.info("Using scoring matrix with affine gap function and following gap penalties: gap_open: {} gap_extend: {}".format(gap_open, gap_extend))
    
    def _sanity_check(self, scoring_matrix, gap_open : int, gap_extend : int):
        lowest_score = min([min(l) for l in scoring_matrix])
        if lowest_score < - gap_open - gap_extend:
            logger.error("Lowest matching score: {} was lower than gap open + gap extend: {}".format(lowest_score, - gap_open - gap_extend))
            raise ValueError

        if gap_open < gap_extend:
            logger.error("Gap extend: {} was greater than gap open: {}.".format(gap_extend, gap_open))
            raise ValueError

    def compare(self, a, b) -> int:
        '''Compares the values a and b and returns this scoring systems match_score if they were equal
        and its mismatch_score if they weren't.

        Parameters:
            a : first value to be compared.
            b : second value to be compared.
        Returns:
            Score resulting from the comparison if the two given values.
        ''' 
        return self.scoring_matrix[self.scoring_alphabet[a]][self.scoring_alphabet[b]]

    def gap_penalty(self, length : int) -> int:
        '''Takes the length of a gap and returns the score penalty for increasing it by one. In this case the score is computed by
        this scoring systems gap_score times the given length.

        Parameters:
            length (int): length of the gap.
        Returns:
            Score computed for increasing the gap with the given length.
        '''
        if length == 1: return -self.gap_open
        elif length > 1: return -self.gap_extend
    