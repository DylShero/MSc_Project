`timescale 1ns / 1ps


module xorshift128 #(
    // default = baseline models default seeds (123456789, 362436069, 521288629, 88675123)
    parameter logic [127:0] SEED_DEFAULT = 128'h075BCD15_159A55E5_1F123BB5_05491333
) (
    input  logic          clk,
    input  logic          rst,       // synchronous, active-high
    input  logic          en,        // advance one word when high
    input  logic          load,      // load seed_in this cycle (overrides en)
    input  logic [127:0]  seed_in,   // starting state {x,y,z,w}, x in MSBs down to w in LSBs (from CPU)
    output logic [31:0]   rnd_out    // current w
);
 
    // One xorshift128 step. from the current state at once, so it is one clock cycle
    //  t = x ^ (x<<11);  new x=y; y=z; z=w;  new w = w ^ (w>>19) ^ (t ^ (t>>8))
    function automatic logic [127:0] xs_step(input logic [127:0] s);
        logic [31:0] x, y, z, w, t, nw;
        {x, y, z, w} = s;
        t  = x ^ (x << 11);
        nw = w ^ (w >> 19) ^ (t ^ (t >> 8));
        return {y, z, w, nw};
    endfunction
 
    // Lock-up guard - the all-zero 128-bit state stays zero forever.
    // Recovers by setting x = 1.
    function automatic logic [127:0] nonzero(input logic [127:0] s);
        return (s == 128'd0) ? {32'h1, 96'h0} : s;
    endfunction
 
    logic [127:0] state;
 
    always_ff @(posedge clk) begin
        if (rst)
            state <= SEED_DEFAULT;             // safe non-zero start
        else if (load)
            state <= nonzero(seed_in);         // CPU starting state
        else if (en)
            state <= xs_step(state);           // next word
    end
 
    assign rnd_out = state[31:0];              // w is the low 32 bits

    
endmodule
