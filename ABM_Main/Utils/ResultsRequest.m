%% Results Request %%

%  This function asks the user if he wants to save the simulation results 
%  stored in the current_area variable. If the answer is 'yes' it stores 
%  the in the ResultsDirectory.

results_answer = input("\n\nSave the Result Matrix? \n 1 = yes \n 2 = no\n\nAnswer: ");

if results_answer == 1
    disp(["Saving ..."]);
    save('Results/current_results.mat', 'tumor_area');
elseif results_answer == 2
    disp(["Ok...\n"]);
else
    disp(["Error ...Close the application!"]);
end