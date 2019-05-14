
--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for schurBetti.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    assert ((schurBetti(3,1,0))#(1,1) == {({4,2}, 1)})
///

TEST ///
    L =  {({30, 30, 15}, 1), ({30, 28, 17}, 1), ({30, 27, 18}, 1), ({30, 26, 19}, 1), 
	  ({29, 26, 20}, 1), ({29, 24, 22}, 1), ({28, 28, 19}, 1), ({28, 27, 20}, 1), 
	  ({28, 26, 21}, 2), ({28, 25, 22}, 1), ({28, 24, 23}, 1), ({27, 26, 22}, 1), 
	  ({27, 24, 24}, 1), ({26, 26, 23}, 1), ({30, 24, 21}, 1), ({30, 25, 20}, 1)};
    assert ((schurBetti(5,2,0))#(13,2) == L)
///

TEST ///
    L = {({14, 4, 4}, 1), ({12, 6, 4}, 1), ({10, 8, 4}, 1), 
	 ({10, 6, 6}, 1), ({8, 8, 6}, 1)};
    assert ((schurBetti(4,2,2))#(5,0) == L)
///


