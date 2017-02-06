<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.reveal {
    color:#06c;
    padding-top: 1em;
    cursor: pointer;
    margin: 5px 0;
}
.answer {
    background-color: #fff;
    padding-left:   1em;
}
</style>

<?php include('faqdata.php'); ?>
<?php
//Output a table of contents with popping up

$sectionCounter = 1;
foreach ($faq as $section) {
  echo '<h4>'.$section['title'].'</h4>';
  echo '<ol>';
    foreach($section['content'] as $item) {
    $label = $item['label'];
  ?>
  <li>
    <div id="<?php echo $label ?>" class="reveal">
      <?php echo $item['question']; ?>
    </div>
    <div id="<?php echo $label ?>_ans" class="answer">
       <?php echo $item['answer'];?>
    </div>
  </li>
  <?php
  }
  echo "</ol><br/>\n";
}
?>
<script>
$(".answer").hide();
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>