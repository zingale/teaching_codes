# compute the odds of getting a given hand in straight 5-card poker.
# M. Zingale
# 2005-03-06

import random

class card:

    def __init__(self, suit=1, rank=2):
        self.suit = suit
        self.rank = rank

    def value(self):
        """ we want things order primarily by rank then suit """
        return self.suit + (self.rank-1)*14

    def __cmp__(self, other):
        return cmp(self.value(), other.value())
    
    def __unicode__(self):
        suits = [u"\u2660",  # spade
                 u"\u2665",  # heart
                 u"\u2666",  # diamond
                 u"\u2663"]  # club
        
        r = `self.rank`
        if self.rank == 11:
            r = "J"
        elif self.rank == 12:
            r = "Q"
        elif self.rank == 13:
            r = "K"
        elif self.rank == 14:
            r = "A"
                
        return r +':'+suits[self.suit-1]

    def __str__(self):
        return unicode(self).encode('utf-8')
        

class deck:
    """ the deck is a collection of cards """

    def __init__(self):

        self.nsuits = 4
        self.nranks = 13
        self.minrank = 2
        self.maxrank = self.minrank + self.nranks - 1

        self.cards = []

        rank = self.minrank
        while (rank <= self.maxrank):

            suit = 1
            while (suit <= self.nsuits):
                currentCard = card(rank=rank, suit=suit)
                self.cards.append(currentCard)

                suit += 1
            rank += 1

    def shuffle(self):
        random.shuffle(self.cards)

    def getCards(self, num=1):
        hand = []

        n = 0
        while n < num:
            hand.append(self.cards.pop())
            n += 1

        return hand
    
    def __str__(self):
        string = ""
        for c in self.cards:
            string += str(c) + " "
        return string


def play():

    nmax = 5000000
    n = 0

    nStraightFlush = 0
    nFourOfAKind = 0
    nFullHouse = 0
    nFlush = 0
    nStraight = 0
    nThreeOfAKind = 0
    nTwoPair = 0
    nOnePair = 0
    
    while (n < nmax):

        # create a deck
        mydeck = deck()
        
        # shuffle
        mydeck.shuffle()
        
        # get a hand
        hand = mydeck.getCards(5)
        hand.sort()
        
        #print hand[0], hand[1], hand[2], hand[3], hand[4]

        found = False
        
        # check for the different hands...

        # straight flush

        # the hand is sorted by rank then suit, make sure
        # that they all have the same suit and that they are
        # sequential
        if (not found and
            (hand[0].suit == \
             hand[1].suit == \
             hand[2].suit == \
             hand[3].suit == \
             hand[4].suit) and
            (hand[0].rank == \
             hand[1].rank - 1 == \
             hand[2].rank - 2 == \
             hand[3].rank - 3 == \
             hand[4].rank - 4)):
            nStraightFlush += 1
            #print "<<< Straight Flush >>>\n"
            found = True


        # four of a kind

        # they are sorted so either cards 0,1,2,3 have the same rank
        # or 1,2,3,4 have the same rank.
        if (not found and
            ((hand[0].rank == hand[1].rank == hand[2].rank == hand[3].rank) or
             (hand[1].rank == hand[2].rank == hand[3].rank == hand[4].rank))):
            nFourOfAKind += 1
            #print "<<< Four of a kind >>>\n"
            found = True

                     
        # full house

        # we are sorted again, so make sure that the first two are equal
        # and then the last three are equal or reverse
        if (not found and
            (((hand[0].rank == hand[1].rank) and
              (hand[2].rank == hand[3].rank == hand[4].rank)) or
             ((hand[0].rank == hand[1].rank == hand[2].rank) and
              (hand[3].rank == hand[4].rank)))):
            nFullHouse += 1
            #print "<<< Full house >>>\n"
            found = True

        
        # flush

        # look for all the same suit
        if (not found and
            (hand[0].suit == \
             hand[1].suit == \
             hand[2].suit == \
             hand[3].suit == \
             hand[4].suit)):
            nFlush += 1
            #print "<<< Flush >>>\n"
            found = True

        
        # straight

        # we are already sorted, so just look at the rank
        if (not found and
            (hand[0].rank == \
             hand[1].rank - 1 == \
             hand[2].rank - 2 == \
             hand[3].rank - 3 == \
             hand[4].rank - 4)):
            nStraight += 1
            #print "<<< Straight >>>\n"
            found = True

            
        # three of a kind

        # since we are sorted, only 0,1,2 or 1,2,3, or 2,3,4 can be
        # equal
        if (not found and
            ((hand[0].rank == hand[1].rank == hand[2].rank) or
             (hand[1].rank == hand[2].rank == hand[3].rank) or
             (hand[2].rank == hand[3].rank == hand[4].rank))):
            nThreeOfAKind += 1
            #print "<<< Three of a kind >>>\n"
            found = True

            
        # two pair and one pair
        if (not found):

            numPairs = 0
            
            if (hand[0].rank == hand[1].rank):
                numPairs += 1

            if (hand[1].rank == hand[2].rank):
                numPairs += 1

            if (hand[2].rank == hand[3].rank):
                numPairs += 1

            if (hand[3].rank == hand[4].rank):
                numPairs += 1

            if numPairs == 2:
                nTwoPair += 1
                #print "<<< Two Pair >>>\n"
                found = True

            elif numPairs == 1:
                nOnePair += 1
                #print "<<< One Pair >>>\n"
                found = True


        
        n += 1


    print "Number of hands: ", nmax
    print " "
    print "  Straight Flush:  ", nStraightFlush, nStraightFlush/float(nmax)
    print "  Four of a kind:  ", nFourOfAKind, nFourOfAKind/float(nmax)
    print "  Full House:      ", nFullHouse, nFullHouse/float(nmax)
    print "  Flush:           ", nFlush, nFlush/float(nmax)
    print "  Straight:        ", nStraight, nStraight/float(nmax)
    print "  Three of a kind: ", nThreeOfAKind, nThreeOfAKind/float(nmax)
    print "  Two pair:        ", nTwoPair, nTwoPair/float(nmax)
    print "  One pair:        ", nOnePair, nOnePair/float(nmax)
    
    
if __name__== "__main__":
    play()


    
        
