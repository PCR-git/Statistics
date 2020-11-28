
# Probability - HW 1
# Peter Racioppo
# 103953689

# (1.)

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random

# Variables
deck_size1 = 52 # Deck size
#sample_size = 10 # Sample Size 1
sample_size = 1000 # Sample Size 2

# This function draws n cards from a deck of size deck_size with replacement
def rand_draw_with_replacement(n, deck_size):
    draw = [] # Initialize array to hold the drawn cards
    i = 1
    # Draw n cards from the deck in the range [0,deck_size-1]
    while i <= n:
      # Select a random card in the range [0,deck_size-1]
      cardi = random.randint(0,deck_size-1)
      draw.append(cardi) # Append card to draw
      i += 1
    #print(draw) # Print array of drawn cards
    return(draw) # Return

# Test: random array of length 10: [20, 40, 42, 13, 32, 49, 7, 4, 3, 50]

# Call rand_draw_with_replacement for n = sample_size and deck_size = deck_size1
vec = rand_draw_with_replacement(sample_size, deck_size1)

# Compute the values "vals" of each of the drawn cards and the number of
# times they were drawn, "cnt"
vals, cnt = np.unique(vec, return_counts=True)
#print(vals) # Print vals
#print(cnt) # Print cnt
sample_freq = cnt/sample_size # Sample frequency
#print(sample_freq) # Print sample frequency

# Plot the sample frequency for each drawn card
fig = plt.figure() # Create figure
ax = plt.axes() # Create plot axes
# Plot histogram of sample frequency
plt.stem(vals,sample_freq,'blue',label='Sample Freq')
# Create array of length deck_size1, with each element of value 1/deck_size1
# This is the probability of selecting any given card
deck_size_vec = [1/deck_size1]*deck_size1
# Plot probability, for comparison
ax.plot(vals,deck_size_vec,'red',label='1/(Deck Size)')
plt.title("Drawing Random Cards with Replacement") # Plot title
plt.xlabel("Card Value") # x-axis label
plt.ylabel("Sample Frequency"); # y-axis label
#plt.legend(); #Legend
plt.show() # Show the plot

# ----------------------------

# (2.)

deck_size1 = 52 # Deck size

# rand_draw2 selects a card from the deck, deletes the card and then selects
# another card from the set of remaining cards.
def rand_draw2(deck_size):
    draw = [] # Initialize array to hold selected cards
    cardi = random.randint(0,deck_size-1) # Draw first card
    draw.append(cardi) # Append to draw
    # Draw second card from a deck with one less card
    cardj = random.randint(0,deck_size-2)
    # In the new deck, every card >= the first selected card becomes indexed
    # by an integer one less than the previous index. We need to map the
    # indexing in deck 2 to the indexing in deck 1. That is, for all
    # cardj >= cardi, we need to increase the indexing by one.
    if cardj >= cardi:
        cardj = cardj + 1
    draw.append(cardj) # Append the second card to draw
    #print(draw) # Print
    return(draw) # Return

# Call rand_draw2 10 times and print the results
i = 1
while i <= 10:
    rand1 = rand_draw2(deck_size1)
    i += 1

# Test: result of calling rand_draw2 10 times
# [11, 1], [1, 8], [31, 44], [44, 8], [14, 28]
# [1, 12], [10, 6], [47, 26], [34, 36], [0, 43]

# rand_draw selects k cards from a deck of size deck_size, without replacement
def rand_draw(deck_size,k):
    # Uses np.random.permutation function to obtain a permuted deck
    vec = np.random.permutation(list(range(deck_size1-1)))
    k_output = vec[0:k] # Selects k first elements of vec
    #print(k_output) # Print
    return(k_output) # Return

k1 = 5; # We will draw 5 permuted cards
# Call rand_draw 10 times and print the results
i = 1
while i <= 10:
    draw_var = rand_draw(deck_size1,k1)
    i += 1

# Test: result of calling rand_draw 10 times
# [11 41 25 35 19]   [19 10 15 43 26]
# [ 4 22 33 37 12]   [32 11 36 40 17]
# [20 36  9 30 50]   [22 41 14 11 43]
# [34 22 47 10 13]   [27 39 42 30 18]
# [11 46 19 17 32]   [34 49 40 38  5]

# ----------------------------

# (3.)

#x1 = [1, 2, 2, 3, 5] # Test vector

# Determines whether there are two pairs in the hand
def is_two_pair(x):
    # Compute the card numbers
    card_num = [] # Initialize array to hold the cards modulo 13

    # Loop over each element in the list and takes its modulus
    for i in range(len(x)):
        card_num.append(x[i] % 13) # Card number modulo 13

    # Returns the unique values and their counts
    vals, cnt = np.unique(card_num, return_counts=True)
    # Returns the values of the counts and the counts of the counts
    vals2, cnt2 = np.unique(cnt, return_counts=True)

# If we interpet the problem to mean that triples count as containing a pair
# and quadruples count as containing two pairs, then we have the following:
# (I'll use this interpretation.)
    if 5 in vals2: # If there are 5 identical elements
        print("Error!") # There should never be 5 of the same element
        num_pair = 0
    elif 4 in vals2: # If there are 4 identical elements, there are 2 pairs
        num_pair = 2
    elif 3 in vals2: # If there are 3 identical elements:
        # If there are 3 identical elements & 2 identical elements, there are 2 pairs
        if 2 in vals2: # There is a triple and a pair
            num_pair = 2
        else: # There is a triple but no pair
            num_pair = 1
    elif 2 in vals2: # If there is at least one pair, the number of pairs is cnt2[1]
        num_pair = cnt2[1]
    else: # Otherwise, there are no pairs
        num_pair = 0

# Alternately, if we interpet the problem to mean that two pairs means that each
# pair has exactly two elements, then we have the following:
#     if 5 in vals2: # If there are 5 identical elements
#         print("Error!") # There should never be 5 of the same element
#         num_pair = 0
#     elif 4 in vals2: # If there are 4 identical elements, there are 0 pairs
#         num_pair = 0
#     elif 3 in vals2: # If there are 3 identical elements:
#         # If there are 3 identical elements & 2 identical elements, there is 1 pair
#         if 2 in vals2: # There is a triple and a pair
#             num_pair = 1
#         else: # There is a triple but no pair
#             num_pair = 0
#     elif 2 in vals2: # If there is at least one pair, the number of pairs is cnt2[1]
#         num_pair = cnt2[1]
#     else: # Otherwise, there are no pairs
#         num_pair = 0

    # true_pair = True if card_num has two pairs
    if num_pair == 2: # There are two pairs
        two_pair = True
    else: # There are less than 2 pairs
        two_pair = False

    #print(two_pair) # Print
    return(two_pair) # Return

# f_pair computes the number of times 2 pairs occur in a draw,
# from a set of N random draws of size k from a deck of size deck_size
def f_pair(deck_size,k,N):
    two_pair_vec = [] # Initialize array to hold True/False values
    i = 1
    while i <= N: # Loop N times
        hand1 = rand_draw(deck_size,k) # Random draw
        pair_parity = is_two_pair(hand1) # Are 2 pairs in the random draw?
        two_pair_vec.append(pair_parity) # True/False value for given draw
        i += 1
    #print(two_pair_vec) # Print
    return(two_pair_vec) # Return

deck_size1 = 52 # Deck size is 52
N1 = 10000 # 10,000 Samples
k1 = 5 # Draw size is 5
truth_vec = f_pair(deck_size1,k1,N1) # Array of N1 true values
# Unique values of truth_vec and their counts
vals, cnt = np.unique(truth_vec, return_counts=True)
# Calculating the number of True and False values, True_num and False_num
if True in vals: # If there are True values
    if False in vals:
        if vals[0] == False: # There are both True and False values
            True_num = cnt[1]
            False_num = cnt[0]
        else: # There are no False values
            True_num = cnt[0]
            False_num = cnt[1]
    else:
        True_num = cnt[0]
        False_num = 0
else: # There are no True values
    True_num = 0
    False_num = cnt[0]

print(True_num/False_num) # Fraction of events in which two pairs are drawn

# We get that the frequency of two-pair hands is about 0.05.

# --------------------------------------
# Theoretical Solution:

# There are (13 choose 2) ways to choose 2 elements from the set of
# 13 elements in the same congruence class.

# There are (4 choose 2) ways to choose 2 members from any set of
# cardinality 4 that is in the same congruence class. Thus, there
# are (4 choose 2)^2 ways to choose 2 members each from 2 sets.

# The remaining element can take on any value other than those already
# selected. That is, it can take on (52-4) values.

# There are 5! ways to order these 5 elements.

# The total number of hands of size 5 drawn from a deck of size 52
# is 52!/(52-5)!

# Thus, the probability of drawing a hand with 2-pairs is:
# [(13 choose 2)*(4 choose 2)^2*5!*(52-4)]/[52!/(52-5)!]
# = 216/4165 ~= 0.052
# Our simulated and theoretical results are in agreement up to
# around the hundreths place.
