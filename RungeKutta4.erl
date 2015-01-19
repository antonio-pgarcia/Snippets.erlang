%% @author antonio.pgarcia
%% @doc @todo Add description to rungekutta.


-module(rungekutta).

%% ====================================================================
%% API functions
%% ====================================================================
-export([main/0, mydydx1/2, mydydx2/2, mydydx3/2, snd/1, rungekutta4/4, rungekutta4/7]).



%% ====================================================================
%% Internal functions
%% ====================================================================

%% --| ODE, y = e ^ x
mydydx1(X, _) ->
	math:exp(X).

%% --| ODE, y = x * e ^ (3 * x) - 2 * y 
mydydx2(X, Y) ->
	X * math:exp(3 * X) - 2 * Y.

%% --| ODE, y = r * y 
mydydx3(_, Y) ->
	R = 0.5,
	R * Y.

%% --| Emulation of Haskell snd (Lisp cdr)
snd(L) ->
	element(2,L).

%% --| Fourth order runge-kutta algorithm calculations
rungekutta4(X, Y, H, F) ->
	K1= H * F(X,Y),
	K2= H * F(X + 0.5 * H, Y + K1 * 0.5),
	K3= H * F(X + 0.5 * H, Y + K2 * 0.5),
	K4= H * F(X + H, Y + K3),
	Y + 1/6 * (K1 + 2 * K2 + 2 * K3 + K4).

%% --| Fourth order runge-kutta solver
rungekutta4(X, Y, H, N, F, XX, YY) ->
	Y1= rungekutta4(X, Y, H, F), 
	if X >= N -> {XX, YY};
	true -> rungekutta4(X+H, Y1, H, N, F, XX ++ [X+H], YY ++ [Y1])
end.

%% --| The main entry point
main() ->
	io:fwrite("Welcome to fourth order Runge-Kutta ODE solver!\n"),
	V1= rungekutta4(0, 1, 0.01, 1.0, fun mydydx1/2, [], []),
	io:fwrite("Value= ~f\n", [lists:last(snd(V1))]),
	V2= rungekutta4(0, 0, 0.01, 1.0, fun mydydx2/2, [], []),
	io:fwrite("Value= ~f\n", [lists:last(snd(V2))]),
	V3= rungekutta4(0, 2, 0.01, 1.0, fun mydydx3/2, [], []),
	io:fwrite("Value= ~f\n", [lists:last(snd(V3))]),
	io:fwrite("end\n").
