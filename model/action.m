function output = action(action_intent, curr_action_focus, xh)

output = action_intent .* curr_action_focus .* xh;

end