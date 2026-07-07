// lfsr32.sv
`timescale 1ns/1ps
module lfsr32 #(
    parameter int unsigned      WIDTH        = 32,
    parameter logic [31:0]      POLY         = 32'h8020_0003,  // taps 32,22,2,1
    parameter logic [31:0]      SEED_DEFAULT = 32'hDEAD_BEEF   // never all-zero
) (
    input  logic              clk,
    input  logic              rst,       // synchronous, active-high
    input  logic              en,        // advance one word when high
    input  logic              load,      // load seed_in this cycle (overrides en)
    input  logic [WIDTH-1:0]  seed_in,   // starting state (from Input from CPU)
    output logic [WIDTH-1:0]  rnd_out    // current LFSR state
);

    function automatic logic [WIDTH-1:0] galois_step(input logic [WIDTH-1:0] s);
        return (s >> 1) ^ ({WIDTH{s[0]}} & POLY);
    endfunction
    
    // Advance a whole word
    function automatic logic [WIDTH-1:0] advance_word(input logic [WIDTH-1:0] s);
        logic [WIDTH-1:0] t;
        t = s;
        for (int i = 0; i < WIDTH; i++)
            t = galois_step(t);
        return t;
    endfunction
    
    // Non-zero guard: all-zero is the LFSR lock-up state (stays zero forever).
    function automatic logic [WIDTH-1:0] nonzero(input logic [WIDTH-1:0] s);
        return (s == {WIDTH{1'b0}}) ? SEED_DEFAULT : s;
    endfunction
    
     logic [WIDTH-1:0] state;
     
     always_ff @(posedge clk) begin
        if (rst)
            state <= SEED_DEFAULT;            // safe non-zero start
        else if (load)
            state <= nonzero(seed_in);        // CPU starting state
        else if (en)
            state <= advance_word(state);     // next word
        // else: hold (flip-flop keeps its value)
    end
    
    assign rnd_out = state; // outputs from always_ff loop to actual logic output
    
endmodule 